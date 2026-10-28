from Controller import *

controller = Controller()

# USING BIN LENGTH OF 1 YEAR (365 DAYS)
# controller.run_all(
#     patient_data_path=r'c:\Users\sarah\RxAlign\Time_Binning_Algorithm\data\test.tsv',
#     med_to_med_class_path=r'c:\Users\sarah\RxAlign\Time_Binning_Algorithm\data\med_class.tsv',
#     bin_length=365,
#     consensus_seq_starts=10,
#     consensus_seq_sub_starts=10,
#     distance_threshold=1,
#     seq_logo_cut_off=0.05,
#     folder=r'c:\Users\sarah\RxAlign\Time_Binning_Algorithm\results\10_patient_example',
# )

# USING TREATMENT LINES
controller.run_all(
    patient_data_path=r'..\data\test.tsv',
    med_to_med_class_path=r'..\data\med_class.tsv',
    use_treatment_lines=True,
    consensus_seq_starts=10,
    consensus_seq_sub_starts=10,
    distance_threshold=0.4,
    seq_logo_cut_off=0.05,
    folder=r'..\results\10_patient_example'
)

# example of query after saving data as json, reloading, and then querying
# query
controller.save_data(folder_path=r'..\results\10_patient_example\instance')

new_controller = Controller()
new_controller.load_data(folder_path=r'..\results\10_patient_example\instance')

query_path = r'..\data\query1.tsv'
new_controller.query(query_file=query_path,
                 med_to_med_class_path=r'..\data\med_class.tsv',
                 dendrogram_path=r'..\results\10_patient_example\query.png',
                 use_treatment_lines=True)