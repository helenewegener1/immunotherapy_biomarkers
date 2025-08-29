# Import the 'logging' module.
import logging as log

# Set the logging options.
log.basicConfig(level = "INFO")

# Import Pandas and the 'ioutil' module from 'bulkDGD'.
import pandas as pd
import bulkDGD
from bulkDGD import ioutil
# Following this: https://bulkdgd.readthedocs.io/en/stable/tutorial_1.html

#######################################################################
################ Step 1 - Preprocess the input samples ################
#######################################################################

date = '2025_08_28'
# dataset = 'Hugo_MEL_counts'
# dataset = 'Riaz_MEL_counts'

datasets = [
    "Auslander_MEL", "Braun_ccRCC", "Cho_NSCLC", "Choueiri_aRCC",
    "Clouphesy_GBM", "Du_MEL", "Gide_MEL", "Hugo_MEL", "Jung_NSCLC",
    "Kim_GC", "Liu_MEL", "Mamdani_rEAC", "Mariathasan_UC",
    "McDermott_RCC", "Miao_ccRCC", "Motzer_aRCC", "Nathanson_MEL",
    "Ravi_NSCLC", "Riaz_MEL", "Rose_mUC", "Snyder_UC",
    "VanAllen_aMEL", "Vandenende_rEAC", "Zappasodi_MEL", "Zhao_GBM"
]

for dataset in datasets:

    # Load the samples into a data frame.
    df_samples = \
    ioutil.load_samples(# The CSV file where the samples are stored
                        csv_file = '/Users/srz223/Documents/projects/project_DGD/cancer_sample_normalization/data/' + dataset + '_counts_dgd_input.csv',
                        # The field separator in the CSV file
                        sep = ",",
                        # Whether to keep the original samples' names/
                        # indexes (if True, they are assumed to be in
                        # the first column of the data frame
                        keep_samples_names = True, 
                        # Whether to split the input data frame into
                        # two data frames, one containing only gene
                        # expression data and the other containing
                        # the extra data about the samples
                        split = False)

    print(df_samples.shape)

    # Preprocess the samples.
    df_preproc, genes_excluded, genes_missing = \
        ioutil.preprocess_samples(df_samples = df_samples)

    print(df_preproc.shape)

    #######################################################################
    ################## Step 2 - Get the trained DGD model #################
    #######################################################################

    # Load the configuration.
    config_model = ioutil.load_config_model("model.yaml")

    # Import the 'model' module from 'bulkDGD.core'.
    from bulkDGD.core import model

    # Get the trained DGD model (Gaussian mixture model and decoder).
    dgd_model = model.BulkDGDModel(**config_model)

    #######################################################################
    ################# Step 3 - Get the optimization scheme ################
    #######################################################################

    # Load the configuration.
    config_rep = ioutil.load_config_rep("two_opt.yaml")

    #######################################################################
    ########### Step 4 - Find and optimize the representations ############
    #######################################################################
    # Get the representations, the corresponding decoder outputs, and
    # the time spent in finding the representations
    df_rep, df_dec_out, df_pred_r_values, df_time_opt = \
        dgd_model.get_representations(\
            # The data frame with the samples - use preprocessed to assure sample has same number of genes as DGD model 
            df_samples = df_preproc,
            # The configuration to find the representations
            config_rep = config_rep)

    df_prob_rep, df_prob_comp = \
        dgd_model.get_probability_density(\
            # The data frame with the samples - use preprocessed to assure sample has same number of genes as DGD model 
            df_rep= df_rep)




    # #######################################################################
    # ###################### Step 5 - Save the outputs ######################
    # #######################################################################

    # Save the preprocessed samples.
    ioutil.save_samples(\
        # The data frame containing the samples
        df = df_preproc,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_preprocessed.csv",
        # The field separator in the output CSV file
        sep = ",")

    # Save the representations.
    ioutil.save_representations(\
        # The data frame containing the representations
        df = df_rep,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_representations.csv",
        # The field separator in the output CSV file
        sep = ",")

    # Save the r-values 
    ioutil.save_representations(\
        # The data frame containing the representations
        df = df_pred_r_values,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_pred_r_values.csv",
        # The field separator in the output CSV file
        sep = ",")

    # Save the decoder outputs.
    ioutil.save_decoder_outputs(\
        # The data frame containing the decoder outputs
        df = df_dec_out,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_decoder_outputs.csv",
        # The field separator in the output CSV file
        sep = ",")

    # Save the time data.
    ioutil.save_time(\
        # The data frame containing the time data
        df = df_time_opt,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_time_opt.csv",
        # The field separator in the output CSV file
        sep = ",")

    ioutil.save_samples(\
        # The data frame containing the samples
        df = df_prob_rep,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_prob_rep.csv",
        # The field separator in the output CSV file
        sep = ",")

    ioutil.save_samples(\
        # The data frame containing the samples
        df = df_prob_comp,
        # The output CSV file
        csv_file = "../out/" + date + "_" + dataset + "_prob_comp.csv",
        # The field separator in the output CSV file
        sep = ",")


