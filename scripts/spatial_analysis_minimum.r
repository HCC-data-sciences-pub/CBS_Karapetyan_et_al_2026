library(readr)
library(dplyr)
library(ggplot2) 
library(spatialTIME)
library(tidyr)


column_id <- "id:mask"
use_example_data <- FALSE


marker_names <- c(
    "classification:status:CD3", "classification:status:HLADR", "classification:status:LAG3", 
    "classification:status:PD1", "classification:status:PDL1", "classification:status:SOX10"
)

################################################################################
################################### Config #####################################

folder_project <- "~/Documents/projects/HCC-CBS-164-Hillman-TBruno"
folder_workflow <- paste(folder_project, "workflows", "workflow_LAG3_2025-03-15_tiles", sep = "/")
folder_spatialtime <- paste(folder_workflow,  "spatialtime", sep ="/")
folder_output <- paste(folder_spatialtime,  "results", sep = "/")
folder_compartments <- paste(folder_spatialtime, "compartments", sep = "/")

get_barcode_from_path <- function(path){
    barcode = strsplit(path, ".", fixed = TRUE)[[1]][1]
    return(barcode)
}
get_barcodes_from_folder <- function(folder){
    filenames <- list.files(folder, pattern = "*.tsv")
    barcodes <- lapply(filenames, get_barcode_from_path)
    return (barcodes)
}

################################################################################
######################### Define Workflow Functions ############################
################################################################################


get_data <- function(folder_project, barcode){
    # Loads the two summary tables and the sample-specific tables for spatialtime
    # Parameters
    #   `folder_project`: str
    #       Path to the project's folder. Only used to load the example data.
    #   `barcode`: str
    #       Barcode of the compartment to load.
    # Returns
    #   List[data.frame, data.frame, List[data.frame]]
    #     - The first element in the list is a dataframe of the clinical data for all samples
    #     - The second element in the list is a dataframe of the summary table
    #     - The third element in the list is a single-item list with the compartment table as a dataframe.
    barcodes <- c(barcode)
    
    sample_filenames <- list()
    
    for (barcode in barcodes){
        name = paste(barcode, "spatial.tsv", sep = ".")
        path = paste(folder_compartments, name, sep ="/")
        sample_filenames <- append(sample_filenames, path)
    }
    
    # These tables should have already been generated when the expression measurement tables were exported to
    # spatialtime's format.
    filename_table_summary  <- paste(folder_workflow, "spatialtime", "summary.tsv",  sep = "/")
    filename_table_clinical <- paste(folder_workflow, "spatialtime", 'clinical.tsv', sep = "/")
    
    # Read in the tables
    table_summary <- read_tsv(filename_table_summary)
    table_clinical <- read_tsv(filename_table_clinical)
    data_spatial <- lapply(sample_filenames, read_tsv)
    data_spatial <- lapply(data_spatial, as.data.frame)
    names(data_spatial) <- c(barcodes)
    output<-list(table_clinical,table_summary, data_spatial)
    return (output)
}


calculate_spatial_interactions <- function(table_clinical, table_summary, data_spatial, marker_names, number_of_permutations = 500, r_range = 0:300){
    # Generates the univariate and bivariate spatialtime tables for a sample
    
    # Convert the input 
    data_spatial <- lapply(data_spatial, as.data.frame)
    xf <-create_mif(
        clinical_data = table_clinical,
        sample_data = table_summary,
        spatial_list = data_spatial,
        patient_id = "id:patient",
        #patient_id = 'id:sample',
        sample_id = 'id:mask'
    )
    
    ################################################################################
    ########################## Univariant Count Method #############################
    ################################################################################
    
    xf <- ripleys_k(mif = xf, mnames =  marker_names,
                    num_permutations = number_of_permutations, method = "K",
                    edge_correction = 'translation', r_range = 0:300,
                    permute = FALSE,
                    keep_permutation_distribution = FALSE, overwrite = TRUE, workers = 4)
    
    ################################################################################
    ########################## Bivariant Count Method ##############################
    ################################################################################
    
    xf <- bi_ripleys_k_WSI(mif = xf, mnames = marker_names,
                           num_permutations = number_of_permutations,
                           permute=FALSE,
                           edge_correction = 'translation', r_range = 0:300,
                           keep_permutation_distribution = FALSE, workers = 4)
    
    ################################################################################
    ##################### Univariant nearest-neighbors Method ######################
    ################################################################################
    xf <- NN_G(mif = xf, mnames = marker_names, num_permutations = number_of_permutations,
               edge_correction = 'rs', r = 0:300, workers = 4)
    
    return(xf)
}

workflow_sample <- function(barcode, folder_output, generate_figures = FALSE){
    path_output_data <- paste(folder_output, paste(barcode, 'xf', 'rds', sep = "."), sep = "/")
    path_output_table_univariant <- paste(folder_output, paste(barcode, 'univariate','tsv', sep = "."), sep = "/")
    path_output_table_bivariant <- paste(folder_output, paste(barcode, 'bivariate','tsv', sep = "."), sep = "/")
    path_output_table_univariant_nn <- paste(folder_output, paste(barcode, 'univariatenn','tsv', sep = "."), sep = "/")
    
    # Load the sample data
    workflow_data <- get_data(folder_project, barcode)
    
    table_clinical <- workflow_data[[1]]
    table_summary <- workflow_data[[2]]
    data_spatial <- workflow_data[[3]]
    
    # Generate the univariate and bivariate tables for a single sample.
    xf <- calculate_spatial_interactions(table_clinical, table_summary, data_spatial, marker_names)
    
    # Save the data
    saveRDS(xf, path_output_data)
    
    write.table(xf$derived$univariate_Count, path_output_table_univariant, sep = "\t", row.names = FALSE)
    write.table(xf$derived$bivariate_Count, path_output_table_bivariant, sep= "\t", row.names = FALSE)
    write.table(xf$derived$univariate_NN, path_output_table_univariant_nn, sep = "\t", row.names = FALSE)
    
    # A few figures just to check the output
    if (generate_figures){
        xf$derived$univariate_Count  %>%
            filter(Marker != 'classification:status:SOX10') %>%
            ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
            geom_line(aes(color = get(column_id)), show.legend = FALSE) +
            facet_wrap(Marker~., scales = 'free') + theme_bw() 
        
        xf$derived$bivariate_Count  %>%
            filter(Anchor == 'classification:status:SOX10',Counted == 'classification:status:CD68') %>%
            ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
            geom_line(aes(color = get(column_id)), show.legend = TRUE) +
            theme_bw()# + scale_color_manual(values = values)
        
        xf$derived$univariate_NN  %>%
            filter(Marker != 'classification:status:SOX10') %>%
            ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
            geom_line(aes(color = get(column_id))) +
            facet_wrap(Marker~., scales = 'free')+ theme_bw() + 
            scale_color_manual(values = c("#1a5fb4", "#26a269", "#e5a50a", "#c64600", "#a51d2d","#613583" ))
    }
    
    return (xf)
}

##################################### load data #########################################
generate_debug_figures <- FALSE
barcodes_masks <- get_barcodes_from_folder(folder_compartments)

for (barcode in barcodes_masks){
    print(barcode)
    xf <- workflow_sample(barcode, folder_output)
}


