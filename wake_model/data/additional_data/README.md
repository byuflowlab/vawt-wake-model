# Additional Data Files for VAWT Wake Model

This folder can be used to include additional data sets to retrain the wake model. Please ensure that the data sets are .csv files and formatted as:

- First row: titles of the data sets (not read into Cross_Validate_EMG.py)
- Following rows: alternating columns of the lateral positions (1) and velocity data sets (2) for each downstream position
- Last rows: make sure that all the column data takes up exactly the same amount of rows; if no more data is reported, use "null" to fill the space

Refer to the data sets in /data/Figshare/VelocityData and /data/Figshare/VorticityData for examples of how data sets should be formatted.
