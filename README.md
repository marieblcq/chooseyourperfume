<p align="center">
  <img src="assets/logo.png" alt="Project Logo" width="400"/>
</p>

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
Choose Your Perfume 
</h1>

<br>


A python🐍 package to help you find your perfect match, ***Choose your perfume*** is a true love story accessible to all❤️.

## 📝 Package description

Developped in 2️⃣0️⃣2️⃣5️⃣ by EPFL students, ***Choose Your perfume*** is a pip-installable python🐍 package which offers perfume recommendations based on the chosen preferences. It regroups a surprising variety of scents🌸 to ensure the highest level of customer satisfaction. Not only can the client choose his prefered scents, he can also use rate his fondness for each of them, in order to receive the most appropriate recommendations based on his wishes. Furthermore, in addition to the main goal⚽ of the package which is giving perfume recommendations, ***Choose Your Perfume*** offers an educational🔍 description of the molecules the client is looking for, including their structures and an indication of how all these structures are similar to each other, in the hopes of developing in the client an interest for chemistry🧪.

## 👩‍💻 Installation

The first step is to create a new environment, to avoid conflict between different packages' dependencies:

```python
conda create -n perfume python=3.10
```

Don't forget to activate the environment!💡

```python
conda activate perfume
```

As this package is pip-installable, chooseyourperfume can be installed like so:

```python
git clone https://github.com/marieblcq/chooseyourperfume.git
cd chooseyourperfume
pip install -e ".[test,doc]"
```


## 🔥 Usage

Now that the environment is activated and the package is installed all that is needed is to run the streamlit app to access the ***Choose Your Perfume***🌸 interface:

```python
streamlit run app_perfume.py
```

If you would like to see a deeper chemical analysis🧪 of the molecules which represent each categories of scent run in the terminal:

```python
python -m src.chooseyourperfume.chem_analysis
```

This should display a colourful graph with each dot representing a molecule and each colour representing a category!

If you want to receive perfume recommendations without using the interface, no problem ! Here is how the package can be used directly in python🐍 code:

```python
from src.chooseyourperfume.logic_cyp import load_data, score_perfumes

perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df = load_data()

#select scents and weights (optional, to be chosen from: 0.1, 0.3, 0.7, 1.0, 1.5, 0.1 being 'its okay' and 1.5 being 'obsessed!')
selected_scents = ['lily', 'lemon', 'mint']
weights = {'lily' : 1.0, 'lemon' : 0.7, 'mint' : 0.3}

results = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights)

print(results.head())
```

Some python🐍 scripts can be condensed into a single line, as you can see here, it is not the case. This code however still manages to show how to quickly leverage the package's main functionality in just a few lines!

Step by step of what this function does💡:
First, it imports from the logic_cyp python🐍 code the main functions needed to obtain a perfume recommendation, which is the purpose of the app
Second step, the data which has been imported needs to be loaded
Third step, the scents🌸 and their weights needs to be selected, this is the part of the code which can be changed to obtain varying perfume recommendations. The scents however have to be chosen from the existing list of scents used in the project.
Final step is to run the code and print the result📝, and hopefully obtain a perfume match which satisfies the demand!

In short, this example demonstrates how to use the core recommendation engine programmatically to get a ranked list of perfumes based on your preferred scent🌸 notes.

## 🛠️ Development installation

If you want to modify this code and add in your own touch, here are a few guidelines to help!

Initialize Git (only for the first time). 

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:mariebilocq/chooseyourperfume.git 
git push -u origin main
```

Then add and commit changes as usual. 

The package should have already been installed by following the instructions above, but here is a reminder just in case (to be run in bash terminal):

```
(perfume) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(perfume) $ python -m pytest --cov=src/chooseyourperfume --cov-report=term --cov-report=html
```

Note: perfume is the example name chosen in this README for the new environment created, it can be modified to your liking.



