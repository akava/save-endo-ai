# Промпт для следующей сессии

Скопируй текст ниже в новую сессию Claude Code (в облаке), в корне этого репозитория.

---

Продолжаем задачу **ICFP Contest 2007 «Morph Endo!»** в репозитории `save-endo-ai` (ветка `main`, remote `git@github.com:akava/save-endo-ai.git`). Нужно найти DNA-префикс с минимальным `risk = 10 × неверные_пиксели + длина_префикса`. Цель — дойти как можно дальше.

## Правила работы (от пользователя)

- Общаемся **по-русски**.
- Решаем **честно**: без интернета, без опубликованных решений и разборов этого контеста. Из своей памяти о контесте тоже не берём готовых префиксов и «знаний»: всё находим экспериментами.
- **Коммитить и пушить после каждого законченного шага**: заработал инструмент, подтвердилась находка, обновлён отчёт. Не копить изменения.
- Вести отчёты в `reports/`:
  - `process.md`: как работали, код, инструменты, техники;
  - `findings.md`: находки по спасению Endo (устройство DNA, флаги, гены, патчи, префиксы, score), **с картинками**;
  - `story.md`: история спасения для девочки 7 лет, с картинками; обновлять при крупных успехах.
- В отчёты класть **только значимые, отобранные картинки** (`reports/img/`). Весь перебор остаётся в `out/` (он в .gitignore).
- Лучшие префиксы сохранять в `prefixes/NN_name.dna` и вести таблицу в findings.md.
- Тяжёлые переборы запускать с таймаутами и в фоне. Пользователь заметил шум вентиляторов, это нормально, но о долгих нагрузках стоит предупреждать.

## Подготовка окружения

```sh
make data/endo.dna && make            # на Linux: make CXX=g++ (нужен zlib: -lz)
python3 -m venv .venv && .venv/bin/pip install pillow numpy
./build/dna -o out/endo.rna           # эталон: 1891886 итераций, 302450 RNA, cost 192646205, ~2 с
./build/build out/endo.rna out/endo.png
```

Сначала прочитай `README.md`, `reports/findings.md` и `reports/process.md`: там всё, что уже известно. Ключевые инструменты:

- `tools/endo.py`: кодировщики (`P()`, `T()`, `nat`, `lit`, `word`), патчи (`set_base`, `write_at`, `no_night`, `call_from_exit`, `stub_gene`), **`combine(*patches)`** для нескольких патчей (каждый патч — функция `off -> str`; skip-адреса должны учитывать длину хвоста префикса!), `run(prefix, name, extra, timeout)`, `score(png)` (приблизительно, по target 300×300).
- `tools/disasm.py DNA start end`: дизассемблер генома (адреса G+n, G = 13615), показывает RNA-команды.
- `tools/rnatree.py RNA depth mindraw`: дерево вызовов генов по RNA (ID-маркеры на входе, `CFPICFP` на выходе).
- `tools/ablate.py RNA GENE_ADDR OUT.png mindraw [iso]`: монтаж абляции или изоляции детей гена.
- `tools/try_genes.py`: вызвать каждый ген после сцены (через exit-ген) и сделать монтаж.
- `analysis/genes.json` (250 генов, [start,end] в G-координатах), `analysis/gene_ids.json` (RNA-ID → ген), `analysis/functions.json`, `analysis/gene_sweep.txt`.
- `build/dna` флаги: `-p/-s` префикс, `-n N`, `-t FROM TO` трасса, `-L file` журнал происхождения (`B` — сравнения, `R` — копирования в координатах исходной DNA), `--dump`.

## Текущее состояние (облачная сессия, 2026-10-02)

Работаем прямо в `main`. Перед каждым коммитом — `CHECKLIST.md` (что обновить: сборка, findings, process, обе истории, картинки, README, NEXT_SESSION). Лучший префикс **`prefixes/55_merged.dna`** (14 124 основания, `tools/build_best.py` собирает его из именованных патчей и адаптеров `tail`, `ecc`, `mu`), **0** неверных пикселей, точно по полноразмерному target `doc/Target-image-600.png` (`score()` в tools/endo.py). История улучшений и таблица — в `reports/findings.md`.

Инструменты: `tools/endo.py` (`kill_instr`, `noop`, `fix_bases`, `push_arg`, `adapter_call`, `crypt_call`, `key128`, `enc_str`), `tools/build_best.py` (`wdiff`, `lit_word_patch`, `retarget_call`, `jmp_at`, `rna_moves`), `tools/polyfit.py` (контур многоугольника по маске, растеризация как у RNA), `tools/hillmodel.py` + `tools/hillridge.py` (точная модель гребней холмов), `tools/polylit.py`, `tools/disasm.py` (`LITMAX`), `tools/flow.py`, `tools/flagrefs.py`, `tools/genetable.py`, `tools/strings.py`, `tools/c/rc4crack2.c`, `tools/grass_emu.py`.

Ключи: «42» (error-correcting-codes), «OPE» (vmu-code), «Out_of_Band_II» (caravan), «9546» (cow-tail), `]` для `goodVibrations` (аудио в `hitWithTheClueStick`), **`no1@Ax3`** (глиф µ, жёлтая записка в гене `sticky`). Ключ травы — исходная фраза плюс `bioMorphPerturb`, который считает биоморф (нужны `enableBioMorph_adaptation` = true и починенный `bioMul`).

Подсказки, которые нашлись в данных: ROT13 на странице безопасности (ключи на жёлтой бумажке), `sticky` (записка), `hitWithTheClueStick` (картинка, PNG, MP3), `shoutOut` (послание пленников: «переставили несколько парабол» — параболы холмов 1 и 2, исправлено), страница Palindromes (зеркальные копии). Эпизоды Major Imp 222 и 285 (чашка под дождём, «смотри на противоположную стену» — чаша кита), `InitialBioMorph.hs` (`enableBioMorph = False`), Beautiful Numbers («the fourth one» → 8128), ImpDoc `checkIntegrity` (сторож для перебора цифрового ключа хвоста).

**Сверка с интернетом** (`reports/audit_search.md`): 9546 — стеганография на портрете E.T. (стр. 112), `sun` = `sunflower` XOR `flower`; поворот 5, контур шарика, текст и холмы другие решатели тоже подбирали по target.

**Аудит чистоты** (`reports/audit_search.md`, проходы 1–6): почти всё, что нашлось перебором, объяснено подсказками или чтением кода; координаты, контур шарика, текст, `mkEmp` 18 и альфа хвоста 7:3 — «чертёж», их честно мерить по target. Открыто: фаза холма 3 (65 → 60) и значение поворота лопастей 5. Проверено и отброшено: неучтённые гены (помощник печати в `printGeneTable`), `sunflower` (испорченная копия `sun`), `transmission-buffer` (послание в духе Arecibo, пасхалка), стеганография, совпадение пятёрок.

## Где остались ошибки (0 px)

1. **Пучки травы** — готово по задуманному пути: `P_biomorph` (`enableBioMorph_adaptation` → `true`, починенный `bioMul`); биоморф пишет `bioMorphPerturb`, который прибавляется к ключу. Обходной `P_grass` (построенное состояние RC4, `tools/rc4craft.py`) больше не нужен.
2. **Фонтан** — готово: сжатая картинка из мёртвого кода `printGeneTable`, зеркальная (словарь cw ↔ ccw), вызывается вместо последнего пучка `grass1` перед китом.
3. **Край воды** — готово (сдвиг только воды: w −= d, s += 2d, фонтан −= 2d).
4. **Длина префикса** (14 124) — теперь весь риск. `tools/build_best.py` собирает через `build_merged` (общая склейка правок). Самые дорогие патчи: шарик ~2600 (литерал многоугольника), облака ~1770, чашка ~1600, `spiro` ~1070, `sun` ~1000. Ключи `crypt` уже копируются из `giveMeAPresent`.
