# save-endo-ai

Решение задачи **ICFP Programming Contest 2007 «Morph Endo!»**, сделанное с нуля по спецификации, без внешних источников.

Задача: найти DNA-префикс, который вместе с DNA пришельца Endo строит картинку, максимально близкую к target:
`risk = 10 × неверные пиксели + длина префикса`.

## Структура

| путь | что там |
|---|---|
| `doc/` | спецификация (`Endo.pdf`), `endo.zip` с DNA, source и target (300×300) |
| `src/dna.cpp` | исполнитель DNA→RNA (персистентный treap, ~2 с на полный прогон Endo) |
| `src/build.cpp` | builder RNA→PNG 600×600 |
| `tools/endo.py` | кодировщики префиксов, запуск, приблизительный score против target |
| `tools/disasm.py` | линейный дизассемблер pattern/template |
| `reports/process.md` | как велась работа: код, инструменты, техники |
| `reports/findings.md` | находки по спасению Endo: устройство DNA, флаги, гены, префиксы |
| `reports/story.md` | история спасения Endo для детей, с картинками |
| `prefixes/` | лучшие найденные префиксы |
| `analysis/` | карта генов (genes.json, gene_ids.json, functions.json), результаты перебора генов |
| `NEXT_SESSION.md` | промпт для продолжения в новой сессии |

## Сборка и запуск

Нужны clang++ (C++17), zlib и Python 3 (pillow и numpy ставятся в `.venv`).

```sh
make data/endo.dna            # распаковать DNA в data/
make                          # build/dna, build/build
python3 -m venv .venv && .venv/bin/pip install pillow numpy

# DNA -> RNA -> PNG
./build/dna -o out/endo.rna                  # без префикса
./build/dna -s IIPIFFCPICICIICPIICIPPPICIIC -o out/hint.rna   # префикс строкой
./build/dna -p prefix.dna -o out/x.rna       # префикс из файла
./build/build out/endo.rna out/endo.png

# префикс + рендер + score одной командой
.venv/bin/python tools/endo.py <PREFIX>
```

### Флаги `build/dna`

| флаг | назначение |
|---|---|
| `-p file` / `-s str` | префикс |
| `-d file` | DNA (по умолчанию `data/endo.dna`; `-` значит пустая) |
| `-o file` | куда писать RNA |
| `-n N` | остановиться после N итераций |
| `-t FROM TO` | трасса итераций в stderr |
| `-L file` | журнал происхождения: проверки (`B`) и копирования (`R`) в координатах исходной DNA |
| `--dump` | вывести оставшуюся DNA в stdout |

### Флаги `build/build`

`build in.rna out.png [--layers prefix] [--stop N]`. `--layers` выгружает все слои с альфой, `--stop N` останавливается после N команд.

## Проверка

Исполнитель воспроизводит эталон из спецификации: 1 891 886 итераций, 302 450 RNA, cost 192 646 205. Префикс из спецификации даёт экран SELFCHECK, на котором все 25 тестов проходят.
