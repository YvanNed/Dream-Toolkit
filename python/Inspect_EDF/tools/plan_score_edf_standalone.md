# Plan : Score_EDF autonome (autoscoring sans le pipeline SomnoLens)

**Statut** : brouillon à relire (2026-10-27). Rien n'est implémenté avant ta confirmation.
**Contexte** : le portage de Score_EDF dans l'outil 3bis est fait (`plan_score_edf_fixes.md`, partie F). Tu veux
aussi garder un outil **autonome** pour qui veut seulement de l'autoscoring, sans passer par les outils 1 à 3 de
SomnoLens (pas de JSON de l'outil 2, pas d'hypnogramme réécrit par l'outil 3).

## Décision à prendre d'abord : corriger Score_EDF, ou le reconstruire sur les bibliothèques de 3bis ?

- **Option A (recommandée) : un Score_EDF « léger » au-dessus des bibliothèques de 3bis.** Le calcul est déjà
  porté, corrigé et vérifié (`autoscoring_lib`, `sleepfm_lib`, `autoscoring_compare_lib` : mêmes résultats que
  Score_EDF à réglages égaux). Le notebook autonome ne garde que ce qui lui est propre : le choix des dossiers,
  la **correspondance des voies par alias** (à la place du JSON de l'outil 2) et la **lecture des hypnogrammes
  bruts** (XML / TXT / CSV, à la place de l'outil 3). Un seul code de calcul à maintenir, et toute correction de 3bis
  profite aux deux.
- **Option B : corriger Score_EDF sur place** (`sleephd/`, `scoredf/`). Deux codes qui divergeront avec le temps.

## Ce que l'outil autonome doit apporter en plus (option A comme B)

1. **Voies par alias, modifiables** : la liste d'alias de Score_EDF étendue aux noms ICM (`EOG G` / `EOG D`,
   `EMG 1` / `EMG 2`, `A1` / `A2`…), dans un petit fichier JSON que l'utilisateur peut compléter, avec le menu de
   correction actuel. Rôles **EOG gauche / EOG droit** (et non `vEOG` / `hEOG`). Détection **par enregistrement**
   (aujourd'hui un choix à la main s'applique à tous les fichiers, même si leur montage diffère).
2. **Hypnogrammes bruts** : XML Compumedics (code 4 = stade 4 → N3, avec avertissement), TXT et CSV avec la table
   de l'outil 3 (`4` → N3, `5` → R, `S1`…`S4`) ; refus clair d'un fichier numérique ambigu au lieu de lire `4`
   comme REM.
3. **Appariement EDF ↔ hypnogramme** : nom exact d'abord, puis préfixe suivi d'un séparateur (aujourd'hui `1.edf`
   peut prendre l'hypnogramme de `117`, et un compagnon `_ScoredEvents_Export.txt` peut être pris pour un
   hypnogramme).
4. **Un modèle qui plante n'arrête plus les autres** sur la même nuit, et chaque modèle n'utilise que les voies
   présentes (corrige le `KeyError: 'hEOG'` et le `KeyError: 'F4'` de ton test).
5. **Mémoire** : lecture des seules voies utiles (corrige le `MemoryError` de la nuit 283 ; pic mesuré à
   16-17 Go sur la nuit 100 avec Score_EDF).
6. **Cache lié aux voies** : changer le mapping ne réutilise plus en silence les anciens résultats.
7. Les corrections déjà faites dans 3bis viennent avec l'option A : référence controlatérale, SleepFM en un seul
   passage + queue courte, recalibrations (division par zéro, N3, threads), `weights_only=False`.

## Bug trouvé pendant le portage (à corriger si l'option B est retenue)

- **Renommage des voies fragile** (`scoredf/channels.py`) : `channel_mapping` remplace temporairement
  `signals._load_edf_channels`, mais `stagers.py` importe cette fonction par son nom. Le renommage ne marche que
  parce que le notebook importe `stagers` *pendant* le premier `with channel_mapping(...)`. Un import fait avant
  (par un autre code, un test) casse tout (`KeyError: 'M1'`, vu dans le script de comparaison).
- Les points vérifiés du rapport (`DreamToolkit_Claude_Report.md`) : #1 (année à 2 chiffres, gravité faible),
  #3, #4, #9, #13, #14 ; #17 et #19 sans objet (voir `plan_score_edf_fixes.md`, phase 2).

## Installation pour un utilisateur sans expérience (partie D du premier plan)

- Un seul installeur `install_autoscore.bat` (crée ou met à jour l'environnement depuis le yml, trouve conda,
  réglage SSL de l'ICM) et un lanceur `autoscore.bat` (remplacent les fichiers `_debugyn`).
- Plus de `SLEEPFM_ROOT` ni de `setup_sleepfm` : le modèle SleepFM se télécharge seul au premier usage
  (`sleepfm_lib.download_sleepfm`, empreintes vérifiées), comme dans 3bis.
- Environnement : celui de 3bis (`autoscore`), aligné si possible sur les versions de SomnoLens (python 3.12,
  numpy 2, torch 2.8 : plus besoin de torch < 2.6), après comparaison de SleepFM dans les deux environnements.

## Questions ouvertes
- Option A ou B ?
- Où vit l'outil autonome : `python/Score_EDF/` (à côté de son code actuel, archivé) ou dans `Inspect_EDF/tools/`
  (même dépôt que les bibliothèques) ?
- Garde-t-on les sorties actuelles de Score_EDF (`autoscore_outputs/`) ou celles de 3bis (par modèle) ?
