# Plan : tool 6, recommander un premier passage brut (A) + enregistrer la préparation du signal (B)

Statut : **à valider**. Le point C (un suffixe de dossier de sortie pour garder un passage brut et un passage
filtré côte à côte) est **reporté** et n'est pas dans ce plan.

## Contexte (2026-10-16)

- En refaisant les figures du poster sur `tools/test_data`, le passage du tool 6 utilisait un passe-haut
  à 0,5 Hz. Rien dans le rapport ne le disait : il a fallu le déduire de la std de C3 (22 µV contre 55 µV brut).
- Mesuré sur les fichiers de défauts, ce passe-haut efface une partie de ce que le tool 6 doit
  détecter : l'écrêtage de 731 n'est plus flaggé (bounds 1,22 % → 0,04 %, hist_extreme 1,21 % → 0 %),
  la dérive de 733 disparaît, la quantification de 738 passe de 3 pics à 1, le segment plat de 735 de 3 pics à 1.
- C'est voulu (SPEC §6, décision du 2026-10-05 : les métriques jugent le signal préparé), et la case est déjà
  décochée par défaut côté EDF. Ce qui manque : (A) dire à l'utilisateur pourquoi faire un premier passage
  brut, et (B) garder la trace de la préparation utilisée.

## Décisions (prises le 2026-10-16)

| # | Question | Décision |
|---|---|---|
| D1 | Format de la provenance dans les TSV | **Une colonne texte `signal_preparation`** : `none` / `high-pass 0.5 Hz` / `resampled to 128 Hz` / `resampled to 128 Hz, high-pass 0.5 Hz`. Case vide = fichier traité avant ce changement (inconnu), distinct de `none` |
| D2 | Skip coché + participant déjà traité avec une autre préparation | **Avertir mais sauter** : le participant reste sauté, une ligne ⚠ nomme les deux préparations |
| D3 | Point C (passages côte à côte) | **Plus tard** |

## A. Texte d'aide sous *Signal preparation* (interface seulement)

- Tool 6 EDF, cellule 3 : un petit texte gris sous le titre *Signal preparation*, avant la ligne d'acquisition :
  > *First run: leave both unticked to see the data as recorded. A high-pass attenuates or erases clipping,
  > quantization, flat segments and DC drift (the metrics are computed on the prepared signal). Use it in a
  > second run, to harmonise a dataset with different acquisition high-passes or to judge the signal as
  > tool 7 will filter it.*
- Twin Curry : le passe-haut y est coché par défaut (0,1 Hz, données DC). `_make_tool6_curry.py` remplace le
  texte par une version Curry, par exemple :
  > *Curry data are DC-coupled, so the high-pass is ticked by default (0.1 Hz) to remove the DC offset. Untick
  > it only to look at the raw offsets. A high-pass attenuates clipping, quantization, flat segments and drift
  > in the metrics.*
- Aucune sortie ne change.

## B. Enregistrer la préparation du signal

### B1. Bibliothèque `tools/quality_lib.py`

- `prepare_signal(...)` renvoie **`(sf, preparation)`** au lieu de `sf`. `preparation` décrit ce qui a
  **réellement** été appliqué : le rééchantillonnage seulement s'il a eu lieu (jamais de suréchantillonnage, et
  une étape qui échoue n'est pas comptée), le passe-haut seulement s'il a réussi. Libellés : `none`,
  `resampled to 128 Hz`, `high-pass 0.5 Hz`, ou les deux séparés par `, `.
- `compute_channel_metrics(..., signal_preparation=None)` : le libellé est écrit dans une nouvelle colonne
  `signal_preparation` de chaque ligne de `{file_id}_quality_metrics.tsv`, **juste après `n_peaks`**
  (`suspect_reason`, `exclude` restent les deux dernières colonnes du calcul). La table par stade ne change pas
  (elle se joint par `file_id`).
- Deux appelants seulement : le tool 6 et l'étape 5 du tool 0. Les deux sont mis à jour. L'affichage du
  tool 0 (« Signal analysed: the EEG as recorded, … ») ne change pas.

### B2. Notebook du tool 6

- **Rapport par participant** : une ligne toujours présente en haut de l'Overview, juste sous le cadre des
  flags, avec la même formulation que le tool 0 :
  *Signal analysed: the EEG as recorded, high-passed at 0.5 Hz (256 Hz), no re-reference.* (ou
  *the EEG as recorded (256 Hz), no re-reference.* sans préparation).
- **`dataset_overview.html`** : une ligne *Signal preparation* dans l'en-tête, qui regroupe les participants par
  libellé, par exemple `none (5) · high-pass 0.5 Hz (3)`. S'il y a plus d'un groupe : un avertissement orange
  (*mixed signal preparations: the metrics are not comparable between these groups*). Les cases vides
  s'affichent comme `not recorded (processed before this was tracked)`.
- **Skip** (D2) : sur la ligne `⤼ … skipped`, on ajoute un ⚠ quand la préparation enregistrée diffère de celle du
  passage courant, par exemple
  `⤼ 734_dead_channel (3/8): skipped (report + data exist) ⚠ made with high-pass 0.5 Hz, this run: none (untick Skip to redo it)`.
  Comparaison sur le **passe-haut** et sur un **rééchantillonnage enregistré**. Limite connue : si le passage
  courant rééchantillonne et que le fichier enregistré ne l'a pas été, on ne peut pas savoir sans relire son
  en-tête s'il était déjà à la fréquence cible : pas d'avertissement dans ce cas (le mélange reste visible dans
  `dataset_overview.html`). Une case vide (fichier ancien) ne déclenche pas d'avertissement par fichier.

### B3. Twin Curry

- Le générateur copie `quality_lib` dans la cellule 1 du twin : la nouvelle colonne et le libellé passent
  tels quels. Ses remplacements de `bounds_pct` visent des lignes que ce plan ne touche pas (vérifié).
- Ajout d'un seul remplacement (le texte d'aide de A), puis on relance `_make_tool6_curry.py`.

## Fichiers touchés

- `tools/quality_lib.py` (B1)
- `tools/6_quality_overview_voila.ipynb` : cellule 3 (A), cellule 4 (rapport + Skip), cellule 1
  (`generate_dataset_overview`)
- `tools/0_live_explore_1file_voila.ipynb` : cellule 10 (nouveau retour de `prepare_signal`, libellé passé à
  `compute_channel_metrics`)
- `tools_curry/_make_tool6_curry.py` + le twin régénéré `tools_curry/6_quality_overview_curry_voila.ipynb`
- `tests/test_variants.py` (test dédié)
- `SPEC.md` (§6, voir plus bas)

## Compatibilité

- Les TSV gagnent une colonne. Tool 7 et tool 0 lisent `quality_summary.tsv` / `_quality_metrics.tsv` par nom de
  colonne : pas d'impact. `generate_dataset_overview` calcule ses statistiques sur une liste explicite de colonnes
  numériques : pas d'impact.
- Goldens : la comparaison accepte les colonnes ajoutées (`snapshot.compare_tables`) et ne compare pas les HTML.
  **Pas de régénération ni d'entrée `INTENDED_CHANGES`** : toutes les valeurs existantes restent identiques.
- Anciens rapports sur disque : leurs `_quality_metrics.tsv` n'ont pas la colonne. Après reconstruction,
  `quality_summary.tsv` a une case vide pour eux (= inconnu), jusqu'à ce qu'ils soient retraités.

## Tests

- Test dédié dans `tests/test_variants.py` (variante rééchantillonnage + passe-haut existante) :
  `signal_preparation == 'resampled to 128 Hz, high-pass 0.5 Hz'` dans `quality_summary.tsv`, la ligne
  *Signal analysed* dans le rapport, puis un second passage sans préparation avec Skip coché, dont la sortie
  doit contenir l'avertissement ⚠.
- Passage par défaut : `signal_preparation == 'none'` partout.
- `test_tool0.py::test_step5_equals_tool6` vérifie déjà que le tool 0 produit exactement les mêmes colonnes
  et valeurs que le tool 6 : il couvre l'étape 5.
- `test_curry_twins.py` : le générateur doit appliquer tous ses remplacements.
- Niveau proposé avant de lancer : **standard** (`-m "not full"`) sur `tool6`, `tool0` et `tool7` (qui lit
  `quality_summary.tsv`). Je te redemanderai au moment de lancer.

## SPEC.md (à valider avec le code)

- §6 *Signal preparation* : le texte d'aide (EDF + Curry), la colonne `signal_preparation` (libellés, case vide =
  inconnu), la ligne *Signal analysed* du rapport, la ligne de `dataset_overview.html`, l'avertissement du Skip
  et sa limite.
- §6 *Outputs per run* : `signal_preparation` dans la liste des colonnes de `quality_summary.tsv`.
- CLAUDE.md : pas de changement (aucune règle transversale nouvelle).

## TODO

- [x] B1 bibliothèque
- [x] B2 notebook du tool 6 (rapport, dataset_overview, Skip)
- [x] A texte d'aide (EDF)
- [x] Tool 0 cellule 10
- [x] B3 générateur Curry + régénération du twin
- [x] Tests dédiés (`test_quality_lib.py` quick, `test_known_answers.py` défaut `none`, variante dans `test_variants.py`)
- [x] Niveau de tests : standard ciblé `(tool0 or tool6 or tool7) and not full` (36 passed) + la variante tool 6
      (passed sur synthetic, skipped sur real comme prévu), 2026-10-18
- [x] SPEC.md + TESTS.md
- [ ] Supprimer ce plan à la fin (commit de suppression séparé)
