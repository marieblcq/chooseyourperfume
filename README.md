![Project Logo](assets/logo.png)

<p align="center">
  <img src="assets/logo.png" alt="Project Logo" width="200"/>
</p>

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
Choose Your Perfume 
</h1>

<br>


A python🐍 package to help you find your perfect match, ***Choose your perfume*** is a true love story accessible to all❤️.

## 📝 Package description

Developped in 2️⃣0️⃣2️⃣5️⃣ by EPFL students, ***Choose Your perfume*** is a pip-installable python🐍 package which offers perfume recommendations based on the chosen preferences. It regroups a surprising variety of scents🌸 to ensure the highest level of customer satisfaction. Not only can the client choose his prefered scents, he can also use rate his fondness for each of them, in order to receive the most appropriate recommendations based on his wishes. Furthermore, in addition to the main goal of the package which is giving perfume recommendations, ***Choose Your Perfume*** offers an educational🔍 description of the molecules the client is looking for, including their names, chemical formulas and lewis drawings, in the hopes of developing in a client an interest for chemistry🧪.

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
pip install .
```

Or via pip and github as well:

```python
pip install git+https://github.com/marieblcq/chooseyourperfume.git
```

## 🔥 Usage

Now that the environment is activated and the package is installed all that is needed is to run the streamlit app to access the ***Choose Your Perfume*** interface:

```python
streamlit run app_perfume.py
```

If you want to receive perfume recommendations without using the interface, no problem ! Here is how the package can be used directly in python🐍 code:

*j'attend que la main function qui donne les recommandations soient finalisée avant de la mettre ici, début de l'exemple*
```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```
The aim of our project is to study the chemistry of perfume!

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.
*fin de l'exemple* 

## 👩‍💻 Installation

*j'ai deja fait ca normalement, je le garde au cas ou, ca commence ici*
Create a new environment, you may also give the environment a different name. 

```
conda create -n chooseyourperfume python=3.10 
```

```
conda activate chooseyourperfume
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(chooseyourperfume) $ pip install jupyterlab
```
*et ca finit la* 

## 🛠️ Development installation

*enfaite de ce que j'ai compris cette partie ne concerne pas les users mais ce sont des indications pour nous, dans les exemples personne n'a gardé cette partie vu qu'elle ne sert pas aux users*
Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:mariebilocq/chooseyourperfume`.

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

To install the package, run

```
(chooseyourperfume) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



