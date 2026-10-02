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

Работаем прямо в `main`. Перед каждым коммитом — `CHECKLIST.md` (что обновить: сборка, findings, process, обе истории, картинки, README, NEXT_SESSION). Лучший префикс **`prefixes/45_parabolas.dna`** (15 357 оснований, `tools/build_best.py` собирает его из именованных патчей и адаптеров `tail`, `ecc`, `mu`), **3 074** неверных пикселя, точно по полноразмерному target `doc/Target-image-600.png` (`score()` в tools/endo.py). История улучшений и таблица — в `reports/findings.md`.

Инструменты: `tools/endo.py` (`kill_instr`, `noop`, `fix_bases`, `push_arg`, `adapter_call`, `crypt_call`, `key128`, `enc_str`), `tools/build_best.py` (`wdiff`, `lit_word_patch`, `retarget_call`, `jmp_at`, `rna_moves`), `tools/polyfit.py` (контур многоугольника по маске, растеризация как у RNA), `tools/hillmodel.py` + `tools/hillridge.py` (точная модель гребней холмов), `tools/polylit.py`, `tools/disasm.py` (`LITMAX`), `tools/flow.py`, `tools/flagrefs.py`, `tools/genetable.py`, `tools/strings.py`, `tools/c/rc4crack2.c`, `tools/grass_emu.py`.

Ключи: «42» (error-correcting-codes), «OPE» (vmu-code), «Out_of_Band_II» (caravan), «9546» (cow-tail), `]` для `goodVibrations` (аудио в `hitWithTheClueStick`), **`no1@Ax3`** (глиф µ, жёлтая записка в гене `sticky`). Не найден: ключ фразы для пучков травы `drawGrassPatch`.

Подсказки, которые нашлись в данных: ROT13 на странице безопасности (ключи на жёлтой бумажке), `sticky` (записка), `hitWithTheClueStick` (картинка, PNG, MP3), `shoutOut` (послание пленников: «переставили несколько парабол» — параболы холмов 1 и 2, исправлено), страница Palindromes (зеркальные копии).

## Где остались ошибки (3 074 px)

1. **Пучки травы** (~2 200 px внизу): target-раскладка известна (`analysis/target_tufts.json`), ключ RC4 не найден. Перебор пользователь запретил: «должна быть подсказка». Важно: в target есть пучки `grass4`, а наш код берёт тип как `b₃ & 3` и `grass4` нарисовать не может; значит, в target отличается и сам код выбора типа, не только ключ.
2. **Фонтан над китом** (~760 px): форма и цвета не найдены ни в одном ведре.
3. **Край воды в чаше** (~120 px).
