# Simulateur Optique Avancé pour Couches Minces

<div align="center">

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://fieldopt.streamlit.app/) &nbsp; [![Python](https://img.shields.io/badge/Python-3.11+-blue?logo=python&logoColor=white)](https://www.python.org/) &nbsp; [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

<br>

<p align="center">
  Une application web développée avec Streamlit pour la simulation et l'optimisation de structures optiques en couches minces.
</p>

<p align="center">
  <img src="https://i.imgur.com/8aZ3v7w.png" alt="Aperçu de l'application"/>
</p>

---

## 📖 À propos du projet

Ce projet est un outil interactif destiné aux ingénieurs, chercheurs et étudiants en optique. Il permet de concevoir, d'analyser et d'optimiser des empilements de couches minces en temps réel, directement depuis un navigateur web.

L'application combine des calculs physiques rapides, basés sur la méthode des matrices de transfert, avec des algorithmes d'optimisation puissants pour aider à la conception de filtres optiques, de miroirs et d'autres composants.

---

## ✨ Fonctionnalités clés

* **Interface Intuitive** : Configurez vos structures complexes en quelques clics.
* **Analyse Spectrale Complète** : Visualisez instantanément la réflectance de votre empilement.
* **Distribution du Champ Électrique** : Analysez en profondeur la répartition de l'énergie lumineuse à l'intérieur de la structure pour une longueur d'onde précise.
* **Optimisation Puissante** :
    * Utilisez des algorithmes d'optimisation globale (`Differential Evolution`) ou locale (`L-BFGS-B`) pour affiner vos designs.
    * Minimisez une fonction coût basée sur la réflectance et des contraintes sur le champ électrique.
* **Export Facile** : Téléchargez vos paramètres de simulation et les résultats spectraux au format `.xlsx`.
* **Calculs Accélérés** : Les routines de calcul les plus intensives sont optimisées avec **Numba** pour une réactivité maximale.

---

## 💻 Technologies utilisées

L'application est construite avec un ensemble d'outils Python modernes et performants :

| Python           | Streamlit        | NumPy            | SciPy            | Matplotlib       | Numba            | Pandas           |
| ---------------- | ---------------- | ---------------- | ---------------- | ---------------- | ---------------- | ---------------- |
| ![Python][py]    | ![Streamlit][st] | ![NumPy][np]     | ![SciPy][sci]    | ![Matplotlib][mp] | ![Numba][nb]     | ![Pandas][pd]    |

[py]: https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white
[st]: https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white
[np]: https://img.shields.io/badge/Numpy-013243?style=for-the-badge&logo=numpy&logoColor=white
[sci]: https://img.shields.io/badge/SciPy-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white
[mp]: https://img.shields.io/badge/Matplotlib-11557c?style=for-the-badge&logo=matplotlib&logoColor=white
[nb]: https://img.shields.io/badge/Numba-00A3E0?style=for-the-badge&logo=numba&logoColor=white
[pd]: https://img.shields.io/badge/Pandas-150458?style=for-the-badge&logo=pandas&logoColor=white

---

## 📈 Suivi d'Utilisation

Pour savoir combien de fois l'application est utilisée, nous pouvons intégrer un service d'analyse simple et respectueux de la vie privée comme [GoatCounter](https://www.goatcounter.com/).

### Étape 1 : Créer un compte GoatCounter

1.  Allez sur [goatcounter.com](https://www.goatcounter.com/) et créez un compte gratuit.
2.  Lors de l'inscription, vous devrez donner un **Code** à votre site. C'est ce code qui sera utilisé dans l'URL (par exemple, `mon-app-optique`). Votre URL GoatCounter sera `https://VOTRE-CODE.goatcounter.com`.

### Étape 2 : Ajouter le script de suivi à Streamlit

Modifiez votre fichier `streamlit_app.py` pour y ajouter le code de suivi. Il est préférable de l'ajouter à la fin du fichier, dans la section de l'interface principale.

```python
import streamlit as st
import streamlit.components.v1 as components

# ... (tout votre code d'application existant) ...

# --- Ajout du compteur de visites ---
st.sidebar.markdown("---")
components.html(
    """
    <script data-goatcounter="[https://VOTRE-CODE.goatcounter.com/count](https://VOTRE-CODE.goatcounter.com/count)" 
            async src="//gc.zgo.at/count.js"></script>
    """,
    height=0,
)
```

**Important** : Remplacez `VOTRE-CODE` par le code que vous avez choisi lors de votre inscription sur GoatCounter.

### Étape 3 : Afficher le compteur sur GitHub (Optionnel)

Pour afficher un badge avec le nombre de visites directement sur cette page, vous pouvez ajouter l'image suivante.

1.  Sur GoatCounter, allez dans *Settings* > *Public dashboard* et activez le tableau de bord public.
2.  Ajoutez l'image suivante à votre `README.md`, en remplaçant `VOTRE-CODE` :

```markdown
<p align="center">
  <img src="[https://VOTRE-CODE.goatcounter.com/counter.svg?path=/](https://VOTRE-CODE.goatcounter.com/counter.svg?path=/)" alt="Compteur de visites"/>
</p>
```

---

## 🚀 Installation et Lancement Local

Pour faire fonctionner ce projet sur votre machine locale, suivez ces étapes :

1.  **Clonez le dépôt**
    ```bash
    git clone [https://github.com/nikonvr/fieldopt.git](https://github.com/nikonvr/fieldopt.git)
    cd fieldopt
    ```

2.  **Créez un environnement virtuel (recommandé)**
    ```bash
    python -m venv venv
    source venv/bin/activate  # Sur Windows: venv\Scripts\activate
    ```

3.  **Installez les dépendances**
    Assurez-vous d'avoir un fichier `requirements.txt` avec le contenu suivant, puis lancez :
    ```bash
    pip install -r requirements.txt
    ```

    Contenu du fichier `requirements.txt` :
    ```text
    streamlit
    numpy
    pandas
    matplotlib
    scipy
    numba
    XlsxWriter
    ```

4.  **Lancez l'application Streamlit**
    ```bash
    streamlit run streamlit_app.py
    ```

L'application devrait maintenant être ouverte dans votre navigateur !

---

## ☁️ Déploiement

Cette application est conçue pour être déployée facilement sur [Streamlit Community Cloud](https://streamlit.io/cloud). Suivez les instructions [ici](https://docs.streamlit.io/streamlit-community-cloud/deploy-your-app) pour déployer votre propre version.
