Simulateur Optique Avancé pour Couches Minces
📖 À propos du projet
Ce projet est un outil interactif destiné aux ingénieurs, chercheurs et étudiants en optique. Il permet de concevoir, d'analyser et d'optimiser des empilements de couches minces en temps réel, directement depuis un navigateur web.

L'application combine des calculs physiques rapides, basés sur la méthode des matrices de transfert, avec des algorithmes d'optimisation puissants pour aider à la conception de filtres optiques, de miroirs et d'autres composants.

✨ Fonctionnalités clés
Interface Intuitive : Configurez vos structures complexes en quelques clics.

Analyse Spectrale Complète : Visualisez instantanément la réflectance de votre empilement.

Distribution du Champ Électrique : Analysez en profondeur la répartition de l'énergie lumineuse à l'intérieur de la structure pour une longueur d'onde précise.

Optimisation Puissante :

Utilisez des algorithmes d'optimisation globale (Differential Evolution) ou locale (L-BFGS-B) pour affiner vos designs.

Minimisez une fonction coût basée sur la réflectance et des contraintes sur le champ électrique.

Export Facile : Téléchargez vos paramètres de simulation et les résultats spectraux au format .xlsx.

Calculs Accélérés : Les routines de calcul les plus intensives sont optimisées avec Numba pour une réactivité maximale.

💻 Technologies utilisées
L'application est construite avec un ensemble d'outils Python modernes et performants :

Python

Streamlit

NumPy

SciPy

Matplotlib

Numba

Pandas

🚀 Installation et Lancement Local
Pour faire fonctionner ce projet sur votre machine locale, suivez ces étapes :

Clonez le dépôt

git clone https://github.com/nikonvr/fieldopt.git
cd fieldopt

Créez un environnement virtuel (recommandé)

python -m venv venv
source venv/bin/activate  # Sur Windows: venv\Scripts\activate

Installez les dépendances
Assurez-vous d'avoir un fichier requirements.txt avec le contenu suivant, puis lancez :

pip install -r requirements.txt

Contenu du fichier requirements.txt :

streamlit
numpy
pandas
matplotlib
scipy
numba
XlsxWriter

Lancez l'application Streamlit

streamlit run streamlit_app.py

L'application devrait maintenant être ouverte dans votre navigateur !

☁️ Déploiement
Cette application est conçue pour être déployée facilement sur Streamlit Community Cloud. Suivez les instructions ici pour déployer votre propre version.
