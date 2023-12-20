version 1.0

workflow sceptre {
    input {
    	String output_directory
        File counts_file
        File design_matrix_file
        File covariates_file
        File grna_targets_file

        #general parameters
        Int cpu = 16
        String memory = "64G"
        String docker = "us.gcr.io/landerlab-atacseq-200218/mehta_sceptre:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call convert_anndata {
        input:
            output_directory = output_directory_stripped,
            counts_file = counts_file,
            design_matrix_file = design_matrix_file,
            covariates_file = covariates_file,
            grna_targets_file = grna_targets_file,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }
    
    output {
        File output_tar = run_sceptre.outputs_tar
    }
}

task run_sceptre {

    input {
        String output_directory
        File counts_file
        File design_matrix_file
        File covariates_file
        File grna_targets_file
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        mkdir -p outputs

        R << CODE
        library(dplyr)
        library(ggplot2)
        library(Matrix)
        library(MatrixExtra)
        library(sceptre)
        library(Seurat)

        # read count matrix and convert to sparse
        counts <- read.csv("~{counts_file}",
                   header=TRUE)
        rownames(counts) <- counts[,1]
        counts <- counts[,colnames(counts) != "X"]
        counts_sparse <- as.csr.matrix(counts, binary = FALSE, logical = FALSE, sort = FALSE)
        counts_sparse <- t(counts_sparse)
        # read design matrix file and convert to sparse
        design_mtx <- read.csv("~{design_matrix_file}",
                   header=TRUE)
        rownames(design_mtx) <- design_mtx[,1]
        design_mtx <- design_mtx[,colnames(design_mtx) != "X"]
        design_mtx_sparse <- as.csr.matrix(design_mtx, binary = FALSE, logical = FALSE, sort = FALSE)
        design_mtx_sparse <- t(design_mtx_sparse)
        # get metadata/covariates file
        covariates <- read.csv("~{covariates_file}",
                   header=TRUE)
        rownames(covariates) <- covariates[,1]
        covariates <- covariates[,colnames(covariates) != "X"]
        # get grna targets and subset to those present in dataset
        grna_targets <- read.csv("~{grna_targets_file}",header=TRUE)
        grna_targets <- grna_targets %>% filter(grna_id %in% rownames(design_mtx_sparse))

        # create sceptre object
        response_matrix <- counts_sparse        # response matrix
        grna_matrix <- design_mtx_sparse        # grna matrix
        extra_covariates <- covariates          # batch information
        grna_target_data_frame <- grna_targets  # gRNA target data frame
        moi <- "low"

        sceptre_object <- import_data(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        grna_target_data_frame = grna_target_data_frame,
        moi = moi,
        extra_covariates = extra_covariates
        )

        sceptre_object

        # get positive controls (on target genes)
        positive_control_pairs <- construct_positive_control_pairs(
        sceptre_object = sceptre_object
        )
        head(positive_control_pairs)

        # construct all sets of discovery pairs
        discovery_pairs <- construct_trans_pairs(
        sceptre_object = sceptre_object,
        positive_control_pairs = positive_control_pairs,
        pairs_to_exclude = "pc_pairs"
        )
        head(discovery_pairs)

        # set analysis parameters
        sceptre_object <- set_analysis_parameters(
        sceptre_object = sceptre_object,
        discovery_pairs = discovery_pairs,
        positive_control_pairs = positive_control_pairs,
        grna_integration_strategy = "singleton",
        side = "left"
        )

        print(sceptre_object)

        # assign guides using maximum method (max umi count)
        sceptre_object <- assign_grnas(
        sceptre_object = sceptre_object, 
        method = "maximum"
        )
        
        #plot outcome of grna assignments
        grnas_to_plot <- c("CD47_1", "CD47_2", "CD47_4")
        pdf(file = "outputs/grna_assignment.pdf", width = 10, height = 10)
        plot(sceptre_object, grnas_to_plot = grnas_to_plot)
        dev.off()

        # qc/plot covariates
        pdf(file = "outputs/covariates.pdf", width = 10, height = 10)
        plot_covariates(sceptre_object)
        dev.off()

        sceptre_object <- run_qc(
        sceptre_object = sceptre_object
        )
        pdf(file = "outputs/qc_plots.pdf", width = 10, height = 10)
        plot(sceptre_object)
        dev.off()

        # calibration check, plot, and results
        sceptre_object <- run_calibration_check(
        sceptre_object = sceptre_object,
        parallel = TRUE
        )
        pdf(file = "outputs/calibration_check.pdf", width = 10, height = 10)
        plot(sceptre_object)
        dev.off()

        calibration_result <- get_result(
        sceptre_object = sceptre_object,
        analysis = "run_calibration_check"
        )
        head(calibration_result)

        # discovery power check
        sceptre_object <- run_power_check(
        sceptre_object = sceptre_object,
        parallel = TRUE
        )
        pdf(file = "outputs/discovery_power.pdf", width = 10, height = 10)
        plot(sceptre_object)
        dev.off()

        # run discovery analysis
        sceptre_object <- run_discovery_analysis(
        sceptre_object = sceptre_object,
        parallel = TRUE
        )
        pdf(file = "outputs/discovery_analysis.pdf", width = 10, height = 10)
        plot(sceptre_object)
        dev.off()

        discovery_result <- get_result(
        sceptre_object = sceptre_object,
        analysis = "run_discovery_analysis"
        )
        head(discovery_result)

        write_outputs_to_directory(
        sceptre_object = sceptre_object, 
        directory = "outputs" 
        )

        CODE

        gsutil -m rsync -r outputs ~{output_dir}
        tar -cvf outputs.tar.gz outputs/
        
    >>>

    output {
        File outputs_tar = "outputs.tar.gz"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(anndata_file, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}