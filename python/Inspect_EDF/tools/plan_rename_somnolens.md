# Plan : renommer `Inspect_EDF` en `SomnoLens` (+ remise à jour des mentions)

Statut : **validé, EN ATTENTE** : à lancer seulement quand tu auras fini ton poster (annotation du 2026-10-16).
Au lancement, je commence par t'envoyer la liste de préparation (Phase 0), tu la fais à la main, puis tu
me dis quand commencer.

## Contexte

- Le dossier `Dream-Toolkit/python/Inspect_EDF/` porte le nom de son premier outil, alors que le pipeline
  couvre maintenant inspection, harmonisation, prétraitement et extraction de mesures, en EDF et Curry.
- Nouveau nom retenu : **SomnoLens**. Vérifié le 2026-10-14 : aucun outil scientifique de ce nom (GitHub,
  PyPI, npm, web). somnolens.com est pris par un site de photographie, et « somnolens » est un mot courant
  (somnolence en danois/suédois/norvégien) : toujours écrire « SomnoLens » + « EEG » ou « sleep ».
- Le dépôt git commence à `Dream-Toolkit/` : le renommage est un `git mv` d'un sous-dossier, l'historique
  est conservé (`git log --follow`).

## Décisions (prises le 2026-10-16)

| # | Question | Décision |
|---|---|---|
| D1 | Renommer l'environnement conda `inspect_edf` | **Oui** → `somnolens` |
| D2 | Casse | **`SomnoLens`** (dossier, titres), **`somnolens`** (environnement, futur paquet) |
| D3 | Fichiers de l'outil 1 (`1_inspect_edf_*`, `1_inspect_curry_voila`) | **Gardés** : le nom décrit bien l'outil |
| D4 | `README.docx` | **Réécrit en `README.md`** (anglais), fonctionnalités actuelles **et futures**, installation légère / complète (voir Phase 6) ; le .docx est supprimé |
| D5 | `README.md` racine de `Dream-Toolkit/` | **Pas réécrit** : on y ajoute seulement un pointeur vers `python/SomnoLens/README.md` |
| D6 | Brouillons non versionnés à chemins absolus (`tools/test_Savitzky-Golay/diagnostics/*.py`) | **Pas touchés** |
| D7 | Fichiers non versionnés (`6bis_metrics_correlation.py`, `inspect_curry.py`, `test_mne_curry_header.py`, `test_data_curry/`, `scientific_com/`) | **Ni versionnés ni documentés** |

## TODO

### Phase 0 : préparation (faite par toi, à mon signal)
- [ ] Je t'envoie cette liste au lancement ; tu me confirmes quand c'est fait.
- [ ] État git propre : tu commites ou mets de côté le travail en cours.
- [ ] Fermer tout ce qui tient des fichiers ouverts dans le dossier : fenêtres VS Code, serveurs Voila/Jupyter,
      noyaux Python, PowerPoint (poster). Sous OneDrive, un fichier verrouillé fait échouer le renommage à moitié.

### Phase 1 : renommer
- [ ] Depuis `Dream-Toolkit/` : `git mv python/Inspect_EDF python/SomnoLens`, puis `git status` doit montrer
      des *renamed* (pas de suppressions + ajouts). Les fichiers non versionnés (test_data, etc.) suivent le dossier.
- [ ] `conda rename -n inspect_edf somnolens` (pas de réseau nécessaire), puis vérifier
      `& "$env:LOCALAPPDATA\miniforge3\envs\somnolens\python.exe" -c "import mne, yasa, gssc"`.
- [ ] Rouvrir VS Code sur `python/SomnoLens/`.

### Phase 2 : mentions actives à mettre à jour (vérifiées par grep le 2026-10-14, à re-greper au lancement)
Dossier (`Inspect_EDF`) :
- [ ] `CLAUDE.md` : titre, « from `Inspect_EDF/` » (section *Canonical launch cwd*).
- [ ] `SPEC.md` : titre, racine de l'arborescence, « repo root (`Inspect_EDF/`) ».
- [ ] `tests/nbdriver.py` (commentaire de `REPO_ROOT`).
- [ ] `tools_curry/_make_tool{2,3,4,6,7}_curry.py` : docstring « Run from Inspect_EDF root ».

Environnement (`inspect_edf` → `somnolens`) :
- [ ] `environment.yml` (`name:`).
- [ ] `CLAUDE.md` (*Running Python*) ; `SPEC.md` (env, procédure notebooks, gssc, *How to run*, tests).
- [ ] `tools/autoscoring_lib.py` : **message affiché à l'utilisateur** quand GSSC manque.
- [ ] Docstrings de lancement : `tests/conftest.py`, `tests/make_golden.py`, `tests/refactor_check.py`,
      `tools/6bis_metrics_correlation.py`, `tools/compare_flat_spectrogram_fix.py`,
      `tools/simple_hypnospectro_yasa_vs_fix.py`, `tools/images/_make_unicorn_gif.py`,
      `tools_curry/_make_tool*_curry.py` (7 fichiers), `tools_curry/test_mne_curry_header.py`.
- [ ] `tools/plan_3bis_autoscoring.md` (plan gardé en référence).

Documents :
- [ ] (D4) `README.md` de SomnoLens (remplace `README.docx`) : présentation, outils actuels (0 à 9, Curry
      expérimental), outils prévus, installation légère / complète, démarrage rapide, structure des sorties.
- [ ] (D5) Un pointeur dans `Dream-Toolkit/README.md`.
- [ ] Description de `SPEC.md` (*Project overview*) : elle dit encore « databases in EDF format » ; élargir à
      inspection, harmonisation, prétraitement, extraction de mesures, EDF + Curry (expérimental).
- [ ] **Toi, à la main** : les tutoriels `.docx`, en particulier `tutorials/1_install_python.docx` (4 mentions
      de `inspect_edf`, la commande périmée `voila inspect_edf_voila.ipynb` → `voila tools/1_inspect_edf_voila.ipynb`,
      **les captures d'écran**, et l'installation légère / complète une fois la Phase 6 tranchée).
- [ ] **Toi, plus tard** : « SomnoLens » sur le poster (tu le modifies directement dans PowerPoint).

Ne pas toucher : `tools/old/`, `tools/jupyter_notebook/` (archives gelées), `.ipynb_checkpoints/`, les sorties
de données (`tools/test_data/summary_inspection/`, `tools/ICEBERG_snippet/`, rapports HTML), qui contiennent
d'anciens chemins absolus (voir Phase 5).

### Phase 3 : mettre l'arborescence de SPEC.md à jour (indépendant du renommage)
- [ ] **Listé mais absent** : `tools/preprocessing_phase1_example_scripts/` → supprimer la ligne.
- [ ] **Versionnés mais absents de l'arborescence** : `README.md` (ex-`README.docx`), `pytest.ini`,
      `tutorials/` (4 guides Word : installation, montage Profusion, export EDF), `tests/outside_chain.py`,
      `tests/refactor_check.py`.
- [ ] Vérifier que les noms de dossiers de sortie cités dans l'arbre et dans *How to run* suivent les
      renommages récents (`preproc_epo/`, `reports_preproc/`, `reports_inspection/`).

### Phase 4 : outils annexes
- [ ] Mes mémoires Claude (4 fichiers citent `Inspect_EDF`/`inspect_edf`) : je les mets à jour. Elles sont
      rangées sous `Dream-Toolkit`, donc le renommage ne les perd pas.
- [ ] **Conversations Claude Code** : elles sont rangées par dossier, dans
      `~/.claude/projects/c--Users-…-Dream-Toolkit-python-Inspect-EDF/` (14 conversations aujourd'hui).
      Pour les retrouver dans `/resume` depuis le nouveau dossier, on les **copie** (sans supprimer l'original)
      dans le dossier équivalent `…-python-SomnoLens/`, que Claude Code crée à la première session lancée
      depuis `SomnoLens`. On teste d'abord avec une seule conversation, puis on copie le reste si elle apparaît.
- [ ] `.claude/settings.local.json` ne contient aucun chemin vers `Inspect_EDF` : rien à faire.

### Phase 5 : vérifications
- [ ] `grep -rn "Inspect_EDF\|inspect_edf"` sur les fichiers actifs : il ne doit rester que les noms de
      fichiers de l'outil 1 et les archives.
- [ ] Tests **quick** (`-m quick`) depuis le nouveau dossier avec le nouvel environnement : le renommage ne
      touche aucun calcul, et `nbdriver` calcule la racine du dépôt à partir de son propre emplacement. Le cache
      de la chaîne (`%TEMP%\dtk_tests`) reste valable (empreintes = contenu des fichiers, pas leurs chemins).
- [ ] Lancer une app Voila à la main (ex. `voila tools/1_inspect_edf_voila.ipynb`).
- [ ] Effet attendu, sans gravité : l'outil 1 compare les fichiers sur leur chemin absolu (`path` de
      `FULL_summary_table_edf.tsv`). Les jeux de données **rangés dans le dépôt** (`tools/test_data`,
      `tools/ICEBERG_snippet`) seront donc ré-inspectés une fois. Les études rangées ailleurs ne changent pas.
- [ ] Proposer le message de commit (tu commites toi-même), puis supprimer ce plan dans un commit à part.

### Phase 6 : réflexions à mener (pas forcément dans ce chantier)
**6a. Prendre une option sur le nom** (possible dès maintenant, indépendamment du renommage) :
- [x] **GitHub** : dépôt `SomnoLens` créé le 2026-10-16 avec un README descriptif et la licence **BSD-3**
      (choix provisoire).
- [ ] ⚠ **Licence à confirmer avec le service juridique / valorisation de l'ICM** : en tant que salarié, le code
      appartient en principe à l'employeur. Demander qui est le titulaire du copyright (ligne `Copyright (c) …`
      du fichier `LICENSE`) et si une licence est imposée (les organismes français recommandent souvent
      **CeCILL-B**, l'équivalent BSD conforme au droit français).
- [ ] **PyPI** : le nom `somnolens` est libre, mais PyPI interdit de réserver un nom avec un paquet vide
      (règle PEP 541 : un nom inutilisé peut être réattribué). On le prend en publiant une **première version
      réelle**, même minimale (par exemple l'inspecteur d'en-têtes EDF), le jour où la distribution est décidée (6b).
- [ ] Vérification de marque (TMview pour l'UE, INPI pour la France) avant une publication officielle.
- [ ] (Optionnel) un nom de domaine : somnolens.com est pris ; .org / .io / .fr à vérifier chez un registraire.

**6b. Comment distribuer les outils** (aujourd'hui, ils vivent dans le dépôt global de l'équipe, ce qui ne
rend pas la distribution très classique ni intuitive). Pistes à comparer, rien n'est décidé :
- [ ] **Dépôt GitHub dédié `SomnoLens`**, extrait du dépôt d'équipe en **gardant l'historique** (`git subtree
      split`) ; le dépôt d'équipe le référence (sous-module, ou simple lien).
- [ ] **Paquet installable** (`pip install somnolens`) avec une commande qui ouvre un menu des outils Voila ;
      les données de test et les archives restent hors du paquet.
- [ ] **Installation Windows « en un clic »** pour les non-programmeurs : un script `.bat` qui crée
      l'environnement et lance l'outil 1, en complément des tutoriels.

**6c. Installation légère / complète** (à trancher avant de réécrire le tutoriel d'installation) :
- [ ] **Légère** = tous les outils sauf l'autoscoring : l'environnement actuel **moins** `pytorch-cpu` et `gssc`.
- [ ] **Complète** = légère **+** `pytorch-cpu` + `gssc` (outil 3bis) **+** ce dont a besoin `python/Score_EDF/`
      (plusieurs Go).
- [ ] ⚠ **Conflit de versions à résoudre** : `Score_EDF/environment_autoscore.yml` demande python 3.10,
      numpy 1.26, yasa ≥ 0.7 et torch < 2.6, alors que SomnoLens tourne sous python 3.12, numpy 2.2 et yasa 0.6.
      Trois options : aligner `Score_EDF` sur SomnoLens (à tester), garder **deux environnements séparés**
      (`somnolens` léger + un environnement d'autoscoring), ou une installation par « extras »
      (`pip install somnolens[autoscoring]`) si la piste du paquet (6b) est retenue.
- [ ] Le message de `autoscoring_lib.py` quand GSSC manque devra pointer vers l'installation complète.
- [ ] ⚠ **Licences à vérifier avant de distribuer l'installation complète** (vérifié sur GitHub le 2026-10-16) :
      | Dépendance | Licence |
      |---|---|
      | MNE-Python, YASA, Voila | BSD-3 (compatibles avec une licence BSD-3) |
      | specparam, edfio | Apache-2.0 (compatibles) |
      | **GSSC** (outil 3bis) | **AGPL-3.0** : copyleft strict, y compris pour un usage à travers un réseau |
      - Le code SomnoLens peut rester en BSD-3 tant que GSSC est une dépendance **installée séparément**, mais
        quiconque distribue SomnoLens **avec** GSSC doit respecter l'AGPL pour l'ensemble.
      - `tools/autoscoring_lib.py` (`gssc_combination_logits`) **reproduit la boucle de scoring de GSSC** : si c'est
        une copie de leur code, ce fichier relève lui-même de l'AGPL. À isoler (module optionnel) ou à
        licencier à part.
      - Vérifier aussi les licences des dépendances de `Score_EDF` (pytorch, `sleephd`, etc.).
      - Faire relire ce point par le service juridique de l'ICM (en même temps que 6a).

## Après validation
Supprimer ce fichier une fois le renommage fait et vérifié (dans un commit à part).
