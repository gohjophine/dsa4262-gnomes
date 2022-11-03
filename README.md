# DSA4262 Project Gnomes

Hi, this is Team Gnomes! This repository is a project done for the module DSA4262 Sense-making Case Analysis Health and Medicine in National University of Singapore.

# Table of Contents

- [Project Description](#Project-Description)
- [Usage](#Usage)
  - [User's Guide](#User's-Guide)
    - [Installation](#Installation)
    - [Software Execution](#Software-Execution)
  - [Developer's Guide](#Developer's-Guide)
    - [Installation](#Installation)
    - [Software Execution](#Software-Execution)
  - [Output](#Output)
- [Contributors](#Contributors)
- [Citation](#Citation)
- [Software License](#Software-License)

# Project Description

Task 1: Develop a machine learning method to identify RNA modifications from direct RNA-Seq data

- Write a computational method that predicts m6A RNA modification from direct RNA-Seq data. The method should be able to train a new new model, and make predictions on unseen test data.

Task 2: Prediction of m6A sites in all SG-NEx direct RNA-Seq samples

- Describe the results and compare them across the different cell lines
- Summarise and visualise your observations.

# Usage

## User's Guide

`make_prediction.py`

- Description: Predicts m6a labels by applying the previously trained Random Forest model (`model.sav`) on unseen test data of direct RNA-seq
- Input:
  - Pathname: Direct RNA-Seq test data
- Output:
  - `prediction.csv`: Prediction of the probability scores for each transcript id, transciption position

Here are all the commands needed that you can copy:

```
git clone https://github.com/gohjophine/dsa4262-gnomes.git
pip install -r requirements.txt
cd dsa4262-gnomes/src
python make_prediction.py ../data/test_data.json ../data/prediction.csv
head ../data/prediction.csv
```

### Installation

1. Clone the git repository onto an AWS machine <br />

```
git clone https://github.com/gohjophine/dsa4262-gnomes.git
```

2. Change your directory to the cloned respository folder <br />

```
cd dsa4262-gnomes
```

3. Install the required python packages <br />

```
pip install -r requirements.txt
```

### Software Execution

1. Change your directory to the `src` folder of the cloned repository <br />

```
cd src
```

2. Run `make_prediction.py` with 2 arguments (path to direct RNA-Seq test data, output path) <br /> `python make_prediction.py <path to direct RNA-Seq test data> <output path>`

```
python make_prediction.py ../data/test_data.json ../data/prediction.csv
```

## Developer's Guide

All Jupyter notebooks can be found in the `src` folder and the data required can be found in the `data` folder.

To follow through the process of the project, you may refer to the notebooks in the following order:

1. `read_json.ipynb`
2. `feature_engineering.ipynb`
3. `resample_normalise.ipynb`
4. `feature_selection.ipynb`
5. `rf_model.ipynb`

The notebooks are combined into two Python files for a more convenient execution, namely:

1. `train_model.py` <br/>

   - Description: Trains a Random Forest model using the direct RNA-Seq data (1st argument) and m6A labels (2nd argument) as its predicted class <br />
   - Input:
     - Pathname: Smaller direct RNA-Seq data <br />
     - Pathname: m6a labels <br />
   - Output:
     - `trial_model.sav`: Trained Random Forest model on smaller dataset

2. `make_prediction.py`
   - Description: Predicts m6a labels by applying the previously trained Random Forest model on unseen test data of direct RNA-seq
   - Input:
     - Pathname: Direct RNA-Seq test data
   - Output:
     - `prediction.csv`: Prediction of each transcript id and transciption position

Here are all the commands needed that you can copy:

```
git clone https://github.com/gohjophine/dsa4262-gnomes.git
pip install -r requirements.txt
cd dsa4262-gnomes/src
python train_model.py ../data/data.json ../data/data/data.info
python make_prediction.py ../data/test_data.json ../data/prediction.csv
head ../data/prediction.csv
```

### Installation

1. Clone the git repository onto an AWS machine <br />

```
git clone https://github.com/gohjophine/dsa4262-gnomes.git
```

2. Change your directory to the cloned respository folder <br />

```
cd dsa4262-gnomes
```

3. Install the required python packages <br />

```
pip install -r requirements.txt
```

### Software Execution

1. Change your directory to the `src` folder of the cloned repository <br />

```
cd src
```

2. Run `train_model.py` with 2 arguments (path to datasets) <br /> `python train_model.py <path to direct RNA-Seq data> <path to m6a labels>`

```
python train_model.py ../data/data.json ../data/data/data.info
```

3. Run `make_prediction.py` with 2 arguments (path to direct RNA-Seq test data, output path) <br /> `python make_prediction.py <path to direct RNA-Seq test data> <output path>`

```
python make_prediction.py ../data/test_data.json ../data/prediction.csv
```

## Output

View the first few lines of the output <br />

```
head ../data/prediction.csv
```

# Contributors

1. Christine Sukardi
2. Goh Jophine
3. Khoo Yueh Leng
4. Kuei Yu Fei
5. Ter Shi Yuen

# Citation

Chen, Ying, et al. "A systematic benchmark of Nanopore long read RNA sequencing for transcript level analysis in human cell lines." bioRxiv (2021). doi: https://doi.org/10.1101/2021.04.21.440736 The SG-NEx data was accessed on 8/9/2022 at registry.opendata.aws/sg-nex-data

The pandas development team. (2022). Pandas-dev/pandas: pandas (Version v1.5.1). Zenodo. https://doi.org/10.5281/zenodo.7223478

Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.

Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.

Lemaître, G., Nogueira, F., & Aridas, C. K. (2017). Imbalanced-learn: A python toolbox to tackle the curse of imbalanced datasets in machine learning. Journal of Machine Learning Research, 18(17), 1–5. Retrieved from http://jmlr.org/papers/v18/16-365.html

Chen, T., & Guestrin, C. (2016). XGBoost: A Scalable Tree Boosting System. In Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (pp. 785–794). New York, NY, USA: ACM. https://doi.org/10.1145/2939672.2939785

Van Rossum, G. (2020). The Python Library Reference, release 3.8.2. Python Software Foundation.

Biggs, J., & Madnani, N. (2022, September 25). Factor Analyzer. program documentation, Educational Testing Service. Retrieved from https://factor-analyzer.readthedocs.io/en/latest/factor_analyzer.html.

# Software License

Copyright 2022 Team Gnomes

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
