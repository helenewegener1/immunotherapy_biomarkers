# Import the 'logging' module.
import logging as log

# Set the logging options.
log.basicConfig(level = "INFO")

studies = ["Auslander_MEL", "Cho_NSCLC", "Choueiri_aRCC",
           "Clouphesy_GBM", "Du_MEL", "Gide_MEL", "Hugo_MEL", "Jung_NSCLC",
           "Liu_MEL", "Mamdani_rEAC", "Mariathasan_UC", #"Kim_GC",
           "McDermott_RCC", "Ravi_NSCLC", "Riaz_MEL", 
           "VanAllen_aMEL", "Vandenende_rEAC", "Zappasodi_MEL", "Zhao_GBM"]


# study = "Kim_GC"

for study in studies:

    # Import the 'ioutil' module from 'bulkDGD'.
    from bulkDGD import ioutil

    #################################################################################################
    ## Step 1 - Get the samples’ representations in latent space and corresponding decoder outputs ##
    #################################################################################################

    # Load the preprocessed samples into a data frame.
    df_samples = \
        ioutil.load_samples(# The CSV file where the samples are stored
                            csv_file = f"../../03_bulkDGD/out/{study}_preprocessed.csv",
                            # The field separator used in the CSV file
                            sep = ",",
                            # Whether to keep the original samples' names/
                            # indexes (if True, they are assumed to be in
                            # the first column of the data frame
                            keep_samples_names = True,
                            # Whether to split the input data frame into
                            # two data frames, one containing only gene
                            # expression data and the other containing
                            # additional information about the samples
                            split = False)

    # Get only the first ten rows.
    # df_samples = df_samples.iloc[:10,:]

    # Load the decoder outputs into a data frame.
    df_dec_out = \
    ioutil.load_decoder_outputs(# The CSV file where the decoder outputs
                                # are stored
                                csv_file = f"../../03_bulkDGD/out/{study}_decoder_outputs.csv",
                                # The field separator used in the CSV
                                # file
                                sep = ",",
                                # Whether to split the input data frame
                                # into two data frame, one containing
                                # only the decoder outputs and the other
                                # containing additional information
                                # about the original samples
                                split = False)

    # Get only the first ten rows.
    # df_dec_out = df_dec_out.iloc[:10,:]

    #################################################################################################
    ############################## Step 2 - Get the trained DGD model ###############################
    #################################################################################################

    # Load the configuration.
    config_model = ioutil.load_config_model("../../03_bulkDGD/script/model.yaml")

    # Import the 'model' module from 'bulkDGD.core'.
    from bulkDGD.core import model

    # Get the trained DGD model (Gaussian mixture model
    # and decoder).
    dgd_model = model.BulkDGDModel(**config_model)

    #################################################################################################
    ####################### Step 3 - Perform differential expression analysis #######################
    #################################################################################################

    # print(dir(dgd_model))

    # Get the r-values.
    r_values = dgd_model._r_values

    # Import the 'dea' module from 'bulkDGD.analysis'.
    from bulkDGD.analysis import dea

    # For each sample
    for sample in df_samples.index:

        # Perform differential expression analysis.
        dea_results, _ = \
            dea.perform_dea(# The observed gene counts for the current
                            # sample
                            obs_counts = df_samples.loc[sample,:],
                            # The predicted means - decoder outputs for
                            # the current sample
                            pred_means = df_dec_out.loc[sample,:],
                            # Which statistics should be computed and
                            # included in the results
                            statistics = \
                                ["p_values", "q_values",
                                "log2_fold_changes"],
                            # The r-values of the negative binomials
                            r_values = r_values,
                            # The resolution for the p-values calculation
                            # (the higher, the more accurate the
                            # calculation; set to 'None' for an exact
                            # calculation)
                            resolution = 1e4,
                            # The family-wise error rate for the
                            # calculation of the q-values
                            alpha = 0.05,
                            # The method used to calculate the q-values
                            method = "fdr_bh")

        # Save the results.
        dea_results.to_csv(# The CSV file where to save the results
                        # for the current sample
                        f"../out/DEA_{study}_{sample}.csv",
                        # The field separator to use in the output
                        # CSV file
                        sep = ",",
                        # Whether to keep the rows' names
                        index = True,
                        # Whether to keep the columns' names
                        header = True)
        


