# Plan : corriger et améliorer Score_EDF (scoring automatique YASA / GSSC / SleepFM)

**Statut** : brouillon à relire et annoter (2026-10-26). Rien n'est implémenté avant ta confirmation.
**Périmètre** : `python/Score_EDF/` (notebook `score_edf_voila.ipynb` + paquets `scoredf/` et `sleephd/`), plus
son installation. L'intégration dans l'outil 3bis de SomnoLens viendra ensuite, dans un plan séparé (phase 5).

**Contexte** : le code de Score_EDF vient d'une étude de **siestes** (`NapWindow`, marqueurs `LIGHTS OFF` /
`REVEIL`, courbes TET). Il a été adapté aux nuits (fenêtre = du premier au dernier epoch scoré), mais plusieurs
hypothèses « sieste » sont restées (voir 2.7). Premier test sur 6 nuits ICM le 2026-10-26 : aucun résultat
(diagnostic en partie B).

---

## Partie A. Introduction de référence

### A1. Les modèles : comment ils ont été entraînés, de quoi ils ont vraiment besoin

| | **YASA 0.7** | **GSSC 0.0.9** | **SleepFM-Clinical** |
|---|---|---|---|
| Référence | Vallat & Walker 2021, *eLife* | Hanna & Flöel 2023, *Front. Neuroinform.* | Thapa et al. 2026, *Nature Medicine* |
| Type | Arbre de décision boosté (LightGBM) sur des caractéristiques calculées à la main | Réseau profond : ResNet 1D (une branche EEG, une branche EOG) + GRU bidirectionnel (contexte de la nuit) | « Modèle de fondation » : transformer pré-entraîné par apprentissage contrastif, puis une tête de stadification (transformer + LSTM bidirectionnel) |
| Données d'entraînement | ~3 000 nuits NSRR (âges, sexes, pathologies variés) | 3 314 nuits NSRR : SHHS2 (936), CFS (730), NCHSDB (665), WSC (983) ; test externe DREEM (80 sujets, consensus de 5 experts) | Pré-entraînement : > 585 000 h de PSG, ~65 000 participants (Stanford Sleep Dataset publié, SHHS exclu) |
| EEG | **1 central obligatoire** référencé à la mastoïde (C4-M1 préféré, C3-M2) | C3, C4, F3, F4 référencés à la mastoïde controlatérale (ou Fpz) ; 1 ou 2 EEG recommandés | De 0 à 10 voies « BAS » (EEG + EOG), noms reconnus dans `channel_groups.json` (517 noms) |
| EOG | Optionnel : EOG gauche (LOC/E1-M2) préféré | **HEOG = EOG gauche − EOG droit** (gauche ou droit seul marche aussi ; VEOG déconseillé) | Fait partie de BAS |
| EMG | Optionnel : menton | Non utilisé | Modalité EMG : 1 à 4 voies |
| Autres | Âge/sexe optionnels (inconnus ici : classifieur sans démographie) | Aucun | Modalités **RESP** (1 à 7 voies : thorax, abdomen, flux, SpO2…) et **EKG** (1 à 2 voies) |
| Prétraitement interne | Sous-échantillonnage à 100 Hz, passe-bande 0,4-30 Hz, z-score robuste sur la nuit. **Ne pas filtrer avant** | Passe-bande 0,3-30 Hz, z-score par voie et par nuit, 85,33 Hz (2 560 points / 30 s) | Passe-bas anti-repliement, interpolation à 128 Hz, z-score par voie et par nuit |
| Résolution de sortie | 30 s | 30 s | **5 s** (Score_EDF agrège en 30 s) |
| Performance publiée | ~85 % d'accord avec le consensus d'experts ; N1 le plus faible | > 80 % d'accord sur tous les jeux de test, y compris une seule voie EEG ou EOG seul | F1 moyen 0,70-0,78 en stadification (tâche secondaire : les auteurs disent qu'il faudrait un modèle spécialisé) |
| Licence | BSD-3 | **AGPL-3.0** | **CC BY-NC 4.0** (pas d'usage commercial) |

Détails utiles :
- **YASA** calcule par epoch de 30 s et par voie : écart-type, IQR, asymétrie, aplatissement, passages par zéro,
  paramètres de Hjorth, puissances absolues et relatives par bande, rapports de puissance, entropie de permutation,
  dimensions fractales. Il les lisse ensuite sur 7,5 min centrées et 2 min passées, puis les normalise sur la nuit
  entière : il lui faut donc **la nuit complète**, pas des morceaux. Fréquence d'échantillonnage > 80 Hz.
- **GSSC** fait un « consensus de permutations » : il passe la nuit avec **chaque combinaison** de voies (chaque EEG
  avec l'EOG, chaque EEG seul, l'EOG seul) et garde, pour chaque epoch, la combinaison la plus sûre. Avec les 4 EEG
  et l'EOG que Score_EDF lui donne, cela fait **9 passages sur la nuit entière**. Le contexte (GRU) utilise toute la nuit.
- **SleepFM** découpe chaque voie en jetons de 5 s à 128 Hz. Une « attention » regroupe les voies d'une même modalité
  (d'où le nombre variable de voies, avec un masque pour les voies absentes). Un transformer travaille ensuite par
  tranches de 5 min, et la tête de stadification (LSTM bidirectionnel) parcourt toute la nuit, jusqu'à 11,4 h.
  Le pipeline officiel exige **au moins une voie reconnue dans chacune des 4 modalités** (BAS, RESP, EKG, EMG).
  Score_EDF ajoute ses propres filtres avant (0,3-35 Hz sur EEG/EOG, 10-100 Hz sur l'EMG) : son commentaire indique
  que sans eux le 50 Hz passe et le modèle ne donne presque plus de probabilité au REM.
- **SleepFM recal A** (recalibration niveau A) : apprend une « température » et un biais pour N1, N2 et REM, sur les
  autres enregistrements (validation croisée en laissant un enregistrement de côté). N3 est fixé à 0 car les siestes
  n'avaient pas de N3. Il faut au moins 2 enregistrements.
- **SleepFM recal B** : réapprend la dernière couche linéaire de SleepFM sur ses caractéristiques figées, avec une
  pénalité choisie par validation croisée imbriquée. La ligne N3 reste figée. Il faut au moins 3 enregistrements (6
  conseillés). Les deux recalibrations utilisent **tes hypnogrammes manuels** : c'est une adaptation supervisée. La
  comparaison reste honnête parce que chaque enregistrement est prédit par des paramètres appris sans lui.

**Ce dont chaque modèle a réellement besoin, et ce que donnent les PSG ICM** (noms vus dans `tools/test_data/100.edf`) :

| Modèle | Minimum | Idéal | Montage ICM |
|---|---|---|---|
| YASA | 1 EEG central-mastoïde | + EOG gauche + EMG menton | C4 + A1 (ou C3 + A2) ; `EOG G` ; EMG menton (`EMG 1` / `EMG 2` ?) |
| GSSC | 1 EEG **ou** 1 EOG | C3/C4/F3 + HEOG (G − D) | C3, C4, F3 + A1/A2 ; `EOG G` − `EOG D` |
| SleepFM | 1 voie par modalité (BAS, RESP, EKG, EMG) | Toutes les voies PSG | Tout est là, **sauf F4** : aujourd'hui Score_EDF exige les 6 dérivations AASM, donc SleepFM est sauté |

### A2. Le notebook : ce qu'il calcule, ce qu'il affiche, comment le lire

**Déroulé**
1. **Section 1** : on choisit le dossier des EDF et le dossier des hypnogrammes (`.edf.XML` Compumedics, `.txt` un
   stade par ligne, ou `.csv`). Les fichiers sont appariés par nom.
2. **Section 2** : la correspondance des voies est détectée sur le premier EDF (modifiable dans les menus), puis on
   coche les modèles à lancer.
3. **Run** : pour chaque enregistrement, lecture de l'hypnogramme ; fenêtre d'analyse = du premier au dernier epoch
   scoré ; lancement des modèles (résultats mis en cache dans `autoscore_outputs/cache/`, à côté du dossier des EDF) ;
   puis les recalibrations.
4. **Tableau des epochs** : epochs de 30 s calés sur le début de la fenêtre. Pour chacun : le stade manuel, les
   probabilités de chaque modèle ramenées à 30 s, et le stade prédit (le plus probable). Sauvé dans
   `epochs_all_stagers.csv`.

**Format d'hypnogramme attendu** (le hypnogramme manuel sert de référence et définit la fenêtre d'analyse)
- **Recherche et appariement** : tous les `.xml`, `.txt` et `.csv` du dossier des hypnogrammes (sous-dossiers
  compris) sont candidats. Un EDF `<id>.edf` prend d'abord `<id>.edf.XML` ou `<id>.xml`, sinon le premier fichier
  dont le nom **commence par** `<id>` (fragile, voir 1.9). Le format est choisi d'après l'extension.
- **Point commun à tous les formats** : un stade par epoch, **le premier epoch commence au début de l'EDF** (aucun
  décalage possible). Les epochs non scorés ou de stade inconnu sont exclus des comparaisons. La fenêtre d'analyse va
  du premier au dernier epoch scoré.
- **`.edf.XML` / `.xml` (export Compumedics Profusion)** : lit `<EpochLength>` (30 s si absent) et la liste
  `<SleepStages>` (un `<SleepStage>` par epoch). Codes : `0` = W, `1` = N1, `2` = N2, `3` = N3, `5` = REM,
  `9` = non scoré. Tout autre code (dont `4`, le stade 4 de R&K) devient « inconnu » (voir 1.8).
- **`.txt`** : un stade par ligne (UTF-8, BOM accepté, majuscules/minuscules indifférentes), epochs de **30 s
  imposés**. Étiquettes reconnues : `W`, `WAKE`, `0` → veille ; `1`, `N1` ; `2`, `N2` ; `3`, `N3` ; `R`, `REM` → REM ;
  et **`4` → REM** (piège, voir 1.8). Toute autre étiquette (`MT`, `?`, stades personnalisés…) = inconnu.
- **`.csv`** : une ligne d'en-tête, et une colonne nommée `stage`, `sleep_stage`, `score`, `hypnogram` ou
  `sleep stage` (casse indifférente) ; mêmes étiquettes que le `.txt`, epochs de 30 s.
- **Compatibilité SomnoLens** : les hypnogrammes réécrits par l'outil 3 (`W`, `N1`, `N2`, `N3`, `R`, + `MT` et les
  stades personnalisés, un par ligne) sont lus correctement (`MT` et les stades personnalisés comptent comme
  inconnus). En revanche, un `.txt` **brut numérique** au format R&K (`0`-`5`) serait mal lu : `4` deviendrait REM au
  lieu de N3, et `5` (REM) deviendrait inconnu. Conseil d'ici là : donner à Score_EDF le `.edf.XML` ou le
  hypnogramme réécrit par l'outil 3.

**Résultats affichés, et comment les lire**
- **3.1 Tableau d'accord** (`agreement_summary.csv`) : pour chaque modèle, le nombre d'epochs, l'accord (% d'epochs
  identiques au manuel), le **kappa de Cohen** (accord corrigé du hasard ; repères : 0,61-0,80 bon, > 0,80
  excellent ; l'accord entre deux experts humains tourne autour de 0,76, soit ~83 % d'epochs identiques), le **F1 par
  stade** (compromis entre « trouver tous les epochs de ce stade » et « ne pas en inventer » ; N1 est toujours le plus
  bas, y compris entre humains) et le kappa par enregistrement.
- **Figure d'accord** (`agreement_by_stager.png`) : à gauche, le F1 par stade et par modèle ; à droite, le kappa de
  chaque enregistrement (un symbole par sujet ; le trait noir = tous enregistrements confondus). Un sujet très en
  dessous des autres pour tous les modèles signale plutôt un problème de signal ou de scoring qu'un problème de modèle.
- **3.2 Matrices de confusion** (`confusion_by_stager.png`) : lignes = stade manuel, colonnes = stade prédit, en % de
  la ligne. La diagonale = bien classé. Confusions typiques : N1 ↔ veille, N1 ↔ N2, N2 ↔ N3.
- **3.3 Grille des hypnogrammes** (`hypnograms_by_stager.png`) : une bande par enregistrement, une ligne pour le manuel
  puis une par modèle, une couleur par stade (gris = pas de donnée). On y voit d'un coup d'œil où les modèles
  divergent du manuel au cours de la nuit.
- **3.4 Hypnodensités** (`hypnodensities_<id>.png`) : pour chaque modèle, la probabilité de chaque stade au cours de
  la nuit (aires empilées), avec l'hypnogramme manuel en bas. Une couleur unique = modèle sûr de lui ; des couleurs
  mélangées = epoch ambigu (souvent les transitions et le N1). Utile pour juger la confiance, pas seulement le stade.

Les couleurs des stades ne suivent pas les conventions de SomnoLens (REM en rouge, etc.) : à harmoniser à l'intégration.

---

## Partie B. Diagnostic du premier test (2026-10-26)

| # | Symptôme | Cause | Correction |
|---|---|---|---|
| B0 | Versions inattendues | Le notebook tournait dans l'env `inspect_edf` (YASA 0.6, sans les dépendances SleepFM) | Lancer par `autoscore_debugyn.bat` |
| B1 | `KeyError: 'hEOG'` sur 6 nuits sur 6 | Les noms ICM (`EOG G`, `EOG D`, `EMG 1`…) ne sont pas dans la liste d'alias ; YASA déclare l'EOG « optionnel » mais le réclame toujours ; et l'erreur de YASA arrête tout l'enregistrement, donc GSSC n'a jamais tourné | 1.2, 1.3, 1.4, 1.5 |
| B2 | `MemoryError` sur la nuit 283 | `_load_edf_channels` lit les **58 voies**, toutes remontées à 512 Hz, puis les copie : ~16 Go pour 9 h | 1.1 |
| B3 | SleepFM serait sauté | Score_EDF exige 17 voies fixes, dont F4, absente du montage ICM | 1.6 |

Contournement en attendant : dans la section 2, mettre `hEOG` → `EOG G` (une seule voie : convient à YASA, et GSSC
accepte un EOG gauche seul) et `vEOG` → `EOG D` (utilisé seulement par SleepFM). Le vrai HEOG (G − D) demande la 1.3.

---

## Partie C. Corrections et améliorations, par phases

### Phase 0. Mesurer avant de toucher
- **0.1 Profilage** sur `tools/test_data/100.edf` (7,85 h, 58 voies) et sur une de tes nuits : temps et mémoire de
  chaque étape (lecture EDF, YASA, GSSC, préparation SleepFM, passage SleepFM principal et « queue », recalibrations).
  C'est la référence de l'audit de vitesse (phase 3).
- **0.2 Référence de non-régression** : sauvegarder les sorties actuelles (CSV du cache, `epochs_all_stagers.csv`)
  sur un enregistrement qui marche (voies mappées à la main), pour vérifier ensuite que les refactorisations sans
  effet attendu donnent les mêmes probabilités.

### Phase 1. Faire marcher Score_EDF sur les PSG ICM (bloquant)
- **1.1 Ne lire que les voies utilisées** (`read_raw_edf(include=…)`), et une seule fois par enregistrement, partagée
  entre les modèles (aujourd'hui : 3 lectures complètes, une par modèle). Corrige B2. *(rapport #9)*
- **1.2 Noms de voies ICM** : ajouter `EOG G` / `EOG D`, les EMG menton, `A1` / `A2` (déjà là). Question Q2 : liste
  en dur, petit fichier JSON modifiable, ou réutiliser dès maintenant le JSON de l'outil 2 (`context_channels`) ?
- **1.3 EOG : rôles « EOG gauche » / « EOG droit »** au lieu de `vEOG` / `hEOG` (dans le code d'origine, `hEOG` est en
  fait une voie droite seule). HEOG dérivé = gauche − droit pour GSSC (comme l'outil 3bis par défaut), EOG gauche
  pour YASA, les deux voies pour SleepFM.
- **1.4 YASA : EOG et EMG vraiment optionnels** : s'ils ne sont pas mappés, YASA tourne sans eux (ils sont déclarés
  optionnels mais le code les réclame).
- **1.5 Isoler les erreurs par modèle** : un `try/except` par modèle et par enregistrement, pour qu'une erreur de
  YASA n'empêche plus GSSC et SleepFM sur la même nuit (règle de gestion d'erreurs du projet).
- **1.6 Exigences de SleepFM = celles du modèle** : au moins une voie par modalité, et les dérivations EEG construites
  avec ce qui existe (ex. F3-M2, C3-M2, C4-M1, O1-M2, O2-M1 sans F4-M1). Le nom du cache doit refléter le montage
  réellement utilisé (aujourd'hui `aasm_all_eeg` en dur).
- **1.7 Cache lié aux voies** : aujourd'hui, changer le mapping réutilise en silence les anciens résultats si
  « Force recompute » n'est pas coché. La clé du cache doit inclure les voies utilisées.
- **1.8 Lecture des stades à harmoniser avec SomnoLens** : dans le XML, le code 4 (stade 4 de R&K) → N3 avec un
  avertissement (aujourd'hui « inconnu ») ; dans les `.txt` / `.csv`, ne plus lire `4` comme REM. Proposition :
  reprendre la table de l'outil 3 (`DEFAULT_MAPPING` : `4` → N3, `5` → R, `S1`…`S4`), et refuser un fichier
  numérique ambigu plutôt que de deviner. Et accepter directement le hypnogramme réécrit par l'outil 3 (choix du
  suffixe, comme les autres outils), ce qui prépare l'intégration à 3bis.
- **1.9 Appariement EDF ↔ hypnogramme** : essayer d'abord les noms exacts (`<id>.edf.XML`, `<id>.xml`), puis un
  préfixe **suivi d'un séparateur**. Aujourd'hui un EDF `1.edf` peut prendre l'hypnogramme de `117`, et un compagnon
  Compumedics (`_ScoredEvents_Export.txt`, `_event_xml.csv`) peut être pris pour un hypnogramme.

### Phase 2. Robustesse (points de ton rapport, vérifiés sur le code)
- **2.1 Monkey-patch de GSSC sans `finally`** (`stagers.py`) : confirmé. Si GSSC plante, sa fonction de vote reste
  remplacée pour la suite. → `try/finally`. *(rapport #4)*
- **2.2 Division par zéro dans la recalibration** (`recalibration.py`, poids de recouvrement) : confirmé. Un epoch
  sans aucun pas de 5 s le recouvrant donne des NaN qui cassent l'optimisation. → exclure ces epochs. *(rapport #3)*
- **2.3 `labels = _`** (`stagers.py`) : confirmé ; risque faible (variable locale) mais illisible. → nom explicite.
  *(rapport #14)*
- **2.4 `PROBABILITY_COLUMNS` défini 5 fois** (stagers, sleepfm, recalibration, individual, + une liste locale dans
  alignment) → une seule constante dans `sleephd/io.py`. *(rapport #13)*
- **2.5 Année à 2 chiffres** (`read_edf_start`) : vérifié, **gravité faible** et non « HIGH ». Python et la norme EDF
  ne diffèrent que pour les années 1969-1984 / 2069-2084, et la date ne sert qu'à afficher des heures. L'exemple du
  rapport (1985 → 2085) est faux : Python donne bien 1985. À aligner quand même sur la lecture d'`edf_io_lib`.
  *(rapport #1)*
- **2.6 `torch.set_num_threads(1)` jamais rétabli** (recal B) : après une recal B, un second « Run » dans la même
  session fait tourner SleepFM et GSSC sur **un seul cœur**. → restaurer la valeur précédente.
- **2.7 Hypothèses « sieste » des recalibrations** : biais N3 fixé à 0 (recal A) et ligne N3 figée (recal B) « car les
  siestes n'ont pas de N3 ». Pour des nuits, libérer N3 quand les enregistrements d'apprentissage en contiennent.
- **2.8 Rien à faire** : figures non fermées *(rapport #19)* : le notebook les ferme déjà après affichage ;
  `individual.py` (XGBoost) *(rapport #17)* n'est pas utilisé par le notebook.

### Phase 3. Vitesse (audit, à chiffrer en phase 0)
- **3.1 Une seule lecture EDF, des seules voies utiles** (= 1.1). C'est aussi un gain de temps : MNE remonte toutes les
  voies lues au rythme le plus rapide (la SpO2 à 1 Hz recalculée à 512 Hz…), d'où le rééchantillonnage FFT vu dans
  ta trace d'erreur.
- **3.2 « Queue » de SleepFM (`cover_tail`)** : pour prédire les dernières minutes (< 5 min) que SleepFM laisse
  tomber, le code relance l'encodeur sur **presque toute la nuit**. Cela double l'étape la plus coûteuse. Pistes :
  ne repasser que les dernières tranches (ex. 30 min), ou ne calculer que la dernière tranche et réutiliser le reste.
  À valider : les prédictions ajoutées doivent rester les mêmes (ou très proches).
- **3.3 Recal B refait tout le passage SleepFM** (`head_features`, sur la nuit et sur la queue) pour récupérer les
  caractéristiques internes. → les capturer pendant le premier passage (un seul passage donne les probabilités et les
  caractéristiques). Gain : un passage complet de plus par nuit.
- **3.4 GSSC : 9 passages par nuit** (consensus de permutations). Option : limiter à 1-2 EEG, comme le recommande
  GSSC. Mesurer l'effet sur l'accord avant d'en faire un réglage par défaut.
- **3.5 Threads et lots** : fixer explicitement le nombre de threads torch ; essayer des lots > 1 pour l'encodeur
  SleepFM (aujourd'hui 1), en surveillant la mémoire.
- **3.6 Parallélisme entre enregistrements** (joblib) : seulement pour YASA (léger), et seulement si la mémoire suit.

### Phase 4. Installation simple + versions récentes
- **4.1 Lever `torch < 2.6`** : `weights_only=False` dans les **4** appels `torch.load` de `sleepfm.py` (poids
  officiels, de confiance) + `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1` avant GSSC (comme l'outil 3bis). Le correctif
  `intel-openmp=2024.1` devient alors inutile (torch 2.8 marche sur ce poste).
- **4.2 Viser python 3.12 / numpy 2 / torch 2.8** (les versions de SomnoLens) : faire tourner SleepFM sur une nuit
  test dans les deux environnements et comparer les probabilités. Si c'est identique (à une tolérance près), aligner le
  yml. GSSC et YASA tournent déjà sur ces versions dans `inspect_edf`. `gssc==0.0.9` reste figé (Score_EDF intercepte
  aussi sa fonction de vote, comme 3bis).
- **4.3 Téléchargement automatique de SleepFM** : détaillé en partie D.
- **4.4 Un installeur unique** `install_autoscore.bat` (crée/met à jour l'env) + le lanceur `autoscore.bat`
  (remplacent les fichiers `_debugyn`).

### Phase 5. Intégration dans l'outil 3bis de SomnoLens (plan séparé, plus tard)
Réutiliser le JSON de l'outil 2 (voies + `context_channels`), l'arborescence d'étude (`derivatives/…`,
`reports_…`), les conventions d'affichage, « Skip already processed », les tests, l'anglais partout. Décider ce
qui migre (SleepFM, recalibrations). La question des licences (AGPL de GSSC, CC BY-NC de SleepFM) est déjà notée dans
`plan_rename_somnolens.md` §6c.

---

## Partie D. Installation pour un utilisateur sans expérience (option « b »)

**Aujourd'hui** (4 étapes) : installer Miniforge → créer l'env (ligne de commande) → `setup_sleepfm_debugyn.bat`
(télécharge SleepFM + variable `SLEEPFM_ROOT`) → lancer.

**Option b : SleepFM se télécharge tout seul au premier usage.** Quand SleepFM est coché et introuvable, le code :
1. télécharge depuis le dépôt officiel, à la version figée (`2bcbae04`), uniquement ce qu'il utilise : 4 fichiers de
   code, 1 fichier de config, 2 fichiers de poids (~28 Mo) ;
2. vérifie l'intégrité de chaque fichier (empreinte SHA-256 écrite dans le code) ;
3. les range dans un dossier cache de l'utilisateur (ex. `%LOCALAPPDATA%\SomnoLens\sleepfm-clinical-2bcbae04\`), en
   affichant la progression dans le notebook ;
4. les fois suivantes, ne retélécharge rien.

`SLEEPFM_ROOT` reste possible pour forcer une copie locale (PC hors ligne). L'outil prévu pour ça, **pooch**, est
celui qu'utilise MNE pour ses jeux de données : il est déjà installé avec `mne-base`. Testé le 2026-10-26 : les
téléchargements depuis GitHub passent sur le réseau ICM. C'est d'ailleurs ce que fait MNE pour ses données ; YASA et
GSSC, eux, livrent leurs poids directement dans leur paquet pip.

- **Avantages** : plus de variable d'environnement, plus de script séparé, une version figée et vérifiée ; et on ne
  redistribue pas les poids SleepFM nous-mêmes (licence CC BY-NC), ils viennent du dépôt officiel.
- **Inconvénient** : il faut internet au premier usage de SleepFM (sinon une copie manuelle).

**Cible : 3 gestes** : installer Miniforge (une fois) → double-clic `install_autoscore.bat` → double-clic
`autoscore.bat`. Plus tard (piste 6b de `plan_rename_somnolens.md`) : un paquet `pip install somnolens[autoscoring]`,
ou un installeur `.exe` qui embarque l'environnement (fonctionne hors ligne, mais pèse 1 à 2 Go).

---

## Partie E. Questions ouvertes
- **Q1 Montage ICM** : comment sont enregistrés `EOG G` / `EOG D` (référencés à A1/A2, ou à autre chose) ? Quel est
  l'EMG menton (`EMG 1`, `EMG 2`, la bipolaire `EMG 1 − EMG 2` ?) ? Tes 6 nuits ont-elles le même montage ?
- **Q2 Noms de voies** : liste en dur étendue, petit JSON modifiable, ou le JSON de l'outil 2 dès maintenant ?
- **Q3 Où travailler** : corriger Score_EDF sur place (phases 1-4), puis l'intégrer à 3bis ? Ou aller directement vers
  3bis ?
- **Q4 Validation** : Score_EDF n'a pas de tests. Proposition : un petit script de comparaison avant/après sur 1-2
  nuits (pas de nouveau cadre de tests), plus ton test dans Voila.
- **Q5 Queue SleepFM (3.2)** : une queue calculée avec un contexte réduit (ex. 30 min) te convient-elle, si les
  prédictions restent proches ?
- **Q6 Sorties** : on garde `autoscore_outputs/` à côté du dossier des EDF jusqu'à l'intégration ?

## Hors périmètre : points du rapport qui concernent d'autres dossiers
- **Inspect_EDF** (à vérifier, puis éventuellement un TODO dans SPEC.md) : #8 `edf_io_lib` ouvre un fichier sans le
  fermer ; #10 `load_events` avale ses erreurs en silence ; #11 `_clock_to_hours` (`macrostructure_lib`) suppose une
  nuit (casse pour une sieste de jour) ; #16 fichiers non fermés dans les outils Curry ; #20 décalages d'en-tête
  difficiles à relire.
- **PSG_scoring / PSG_pipeline** (#2, #5, #6, #7, #12, #15, #18) : anciens scripts, non touchés.

## Après validation
Supprimer ce fichier une fois le travail fait et vérifié (dans un commit à part).
