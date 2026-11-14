# AGENT HANDOVER — SPHEREx SkySurvey Backend (`skysurvey-backend`)

Этот документ служит точкой входа и передачи контекста между AI-агентами, разрабатывающими бэкенд астрофизической платформы поиска коричневых карликов, быстродвижущихся объектов и гипотетической Девятой планеты («Планеты X») по данным космического телескопа **NASA SPHEREx**.

---

## 1. Текущее состояние проекта и Git-ветки

- **Репозиторий:** `d:\projectsgit\hakaton\back_spherex\skysurvey-backend`
- **Активная ветка разработки:** `feat/roadmap-phase2` (ответвлена от `main` на коммите `07dcba1`)
- **Статус Phase 2 (Задачи 1–5):** полностью реализованы, покрыты unit/integration-тестами и оформлены отдельными атомарными коммитами:
  - `081d7c4` — `chore(repo): add .gitignore and untrack bytecode caches`
  - `5207870` — `feat(wise): add CatWISE2020/AllWISE IR cross-validation and 10yr proper motion`
  - `c5777f1` — `feat(solar-system): add IMCCE SkyBoT minor planet filtering per epoch`
  - `413cbd0` — `feat(spectrum): add automated Planck Teff fitting and CH4/H2O absorption depth`
  - `fa0de2a` — `feat(astrometry): add multi-epoch parallax vs linear proper motion fitter`
  - `8853c07` — `feat(deploy): add warmup_presets script and Render/Railway production configs`
- **Документация синхронизирована:** `docs/FRONTEND_GUIDE.md`, `docs/AGENT_HANDOVER.md` и корневой `AGENT_HANDOVER.md` (в `d:\projectsgit\hakaton\AGENT_HANDOVER.md`, вне Git-репозитория) поддерживаются идентичными.

### Правило работы для следующих агентов
1. Всегда выполняйте новые изменения в **отдельной ветке** (например, от `feat/roadmap-phase2` или после её слияния в `main`), выделяя **каждую задачу в отдельный атомарный коммит**.
2. После завершения каждой задачи обновляйте `docs/AGENT_HANDOVER.md` и копируйте его в `d:\projectsgit\hakaton\AGENT_HANDOVER.md`, а также обновляйте `docs/FRONTEND_GUIDE.md` при изменении API-схем.

---

## 2. Быстрый старт и архитектура

### 2.1. Локальное окружение и запуск тестов (Windows)
В репозитории настроено виртуальное окружение `.venv` (Python 3.14).

```powershell
cd d:\projectsgit\hakaton\back_spherex\skysurvey-backend

# Запуск сервера
.\.venv\Scripts\uvicorn.exe app.main:app --host 0.0.0.0 --port 8000 --reload
```

> [!TIP]
> **Важно для Windows / песочницы при запуске pytest:**
> В `.env` задано `CACHE_DIR=/tmp/fits_cache`, что на Windows без Docker может приводить к `PermissionError` при инициализации `app.main`. Кроме того, стандартный системный `temp` может вызывать `PermissionError` внутри `asyncio.to_thread`. Запускайте тесты с явным локальным `CACHE_DIR` и `--basetemp` (папка `.pytest_cache/` уже добавлена в `.gitignore`):
>
> ```powershell
> $env:CACHE_DIR="d:\projectsgit\hakaton\back_spherex\skysurvey-backend\.pytest_cache\fits_cache"
> .\.venv\Scripts\python.exe -m pytest tests\test_pipeline.py -v --basetemp="d:\projectsgit\hakaton\back_spherex\skysurvey-backend\.pytest_cache\tmp"
> ```

### 2.2. Реализованные модули (v1.0.0 + Phase 2)

```text
skysurvey-backend/
├── .env / .env.example
├── .gitignore
├── Dockerfile                  # python:3.12-slim, копирует app/, sql/, scripts/, порт ${PORT:-8000}
├── railway.json                # Blueprint конфигурации деплоя Railway
├── render.yaml                 # Blueprint конфигурации деплоя Render
├── pyproject.toml
├── requirements.txt
├── app/
│   ├── config.py               # Settings (IRSA, Gaia, WISE, SkyBoT, Supabase, Cache, Pipeline)
│   ├── main.py                 # FastAPI create_app(), lifespan, /health и роуты v1
│   ├── api/v1/targets.py       # /target, /target/spectrum, /target/diff, /target/restretch, /presets
│   ├── core/
│   │   ├── keys.py             # Детерминированные ключи кеша
│   │   └── spherex.py          # Константы SPHEREx (6 детекторов D1..D6, 102 канала LVF, DEFAULT_FLAG_BITS)
│   ├── schemas/
│   │   ├── spectrum.py         # SpectralModelPoint, SpectrumSummary, ChartPoint, SpectrumResponse
│   │   └── target.py           # TargetResponse, CrossMatch, WiseSource, EpochOut, MinorPlanetOverlay, DiffAnomaly, MotionModel
│   └── services/
│       ├── astrometry_fitter.py    # Task 4: EpochCentroid, measure_subpixel_centroid, fit_motion_model (BIC: линейное PM vs параллакс)
│       ├── cache_manager.py        # FitsCache (LRU + TTL локальный кеш FITS-фрагментов)
│       ├── fits_processor.py       # DisplayGrid, coadd_epoch, deep_stack, compute_difference_map, detect_diff_anomalies
│       ├── fits_stream.py          # Потоковое чтение/распаковка SPHEREx FITS без полной загрузки файла
│       ├── gaia_service.py         # Gaia DR3 TAP + перенос координат с собственным движением (ang_sep_arcsec)
│       ├── irsa_service.py         # Поиск и метаданные SPHEREx L2 Cutouts через IRSA TAP / SIA
│       ├── pipeline.py             # PipelineService — оркестратор этапов, кеша Supabase, кросс-матча и фона
│       ├── solar_system_service.py # Task 2: SolarSystemService (IMCCE SkyBoT на mjd_mid каждой эпохи)
│       ├── spectrum_extractor.py   # Task 3: extract_spectrum, fit_spectral_model (Planck Teff + CH4/H2O)
│       ├── supabase_service.py     # Supabase DB + Storage (с локальным файловым фолбэком)
│       └── wise_service.py         # Task 1: WiseService (CatWISE2020 / AllWISE + propagate_wise 10-yr PM)
├── docs/
│   ├── AGENT_HANDOVER.md
│   └── FRONTEND_GUIDE.md
├── scripts/
│   ├── check_supabase.py       # Проверка подключения к Supabase
│   └── warmup_presets.py       # Task 5: Прогрев PRESET_TARGETS (wise_0855, barnards_star, eps_indi_b, rw_aurigae, r_draconis)
├── sql/001_schema.sql
└── tests/test_pipeline.py
```

#### Краткое описание реализованных фич Phase 2 (Задачи 1–5)
1. **CatWISE2020 / AllWISE Cross-Validation (`app/services/wise_service.py`):**
   - `WiseService.query` запрашивает ИК-каталог `CatWISE2020` через IRSA TAP (с автоматическим фолбэком на `AllWISE`), учитывает, что `catwise_2020` хранит `pmra`/`pmdec` в `arcsec/yr` (конвертация в `mas/yr`), и через `propagate_wise` переносит координаты от эпохи WISE (`CATWISE_EPOCH_MJD = 57170.0` / `ALLWISE_EPOCH_MJD = 55400.0`) на `mjd` наблюдения SPHEREx.
   - `WiseService.best_match` формирует `WiseSource` (`w1mpro`, `w2mpro`, `w1_w2_color`, `is_ultracool_color` при `W1 - W2 >= 1.5`, `baseline_pm_arcsec_yr`), записываемый в `CrossMatch.wise_match`.
2. **IMCCE SkyBoT Minor Planet Filtering (`app/services/solar_system_service.py`):**
   - `SolarSystemService.query_epoch` и `query_epochs` параллельно запрашивают IMCCE SkyBoT на `mjd_mid` каждой эпохи (поддержка JSON и pipe-delimited текста, а также десятичных и sexagesimal координат), проецируют координаты на `512x512` PNG через `DisplayGrid.world_to_png` и заполняют `EpochOut.known_minor_planets` (`MinorPlanetOverlay`).
3. **Planck $T_{\text{eff}}$ Fitting & $\text{CH}_4/\text{H}_2\text{O}$ Absorption Depths (`app/services/spectrum_extractor.py`):**
   - `fit_spectral_model` подбирает масштабированную функцию Планка $B_\nu(\lambda, T_{\text{eff}})$ в диапазоне `[200, 6000]` K вне полос поглощения, оценивает глубины полос `ch4_absorption_depth` (`3.2–3.4 мкм`) и `h2o_absorption_depth` (`2.55–2.85 мкм`) относительно именно кривой Планка (чтобы крутой наклон холодного тела не давал ложного поглощения при линейной интерполяции) и строит сглаженную кривую `model_curve` (`list[SpectralModelPoint]`).
4. **Multi-Epoch Parallax vs Linear Proper Motion Fitter (`app/services/astrometry_fitter.py`):**
   - Для `moving_dipole` при наличии `>= 3` эпох `measure_anomaly_centroids` и `measure_subpixel_centroid` измеряют субпиксельные центроиды по всем доступным коаддам.
   - `fit_motion_model` вычисляет факторы годичного параллакса (`_parallax_factors` через `astropy.coordinates.get_sun`) и сравнивает по критерию BIC линейное собственное движение (4 параметра), чистый параллакс (3 параметра) и при `N >= 4` комбинированный дрейф + параллакс (5 параметров), заполняя `DiffAnomaly.motion_model` (`"linear_proper_motion" | "parallax_dominated" | "ambiguous"`) и `DiffAnomaly.parallax_arcsec`.
5. **Preset Warmup & Production Deploy (`scripts/warmup_presets.py`, `Dockerfile`, `render.yaml`, `railway.json`):**
   - `scripts/warmup_presets.py` экспортирует `warmup_all(...)` и CLI (`--force`, `--target`, `--dry-run`) для синхронного прогрева всех `PRESET_TARGETS` (`wise_0855`, `barnards_star`, `eps_indi_b`, `rw_aurigae`, `r_draconis`).
   - `Dockerfile` обновлён до `python:3.12-slim` с динамическим `${PORT:-8000}` и копированием `scripts/`; добавлены манифесты `render.yaml` и `railway.json` с healthcheck `/health`.

---

## 3. Дорожная карта для следующих агентов (Phase 3)

Ниже перечислены 4 приоритетные научные и продуктовые задачи для развития платформы.

---

### Задача 6 (Высокий приоритет): PSF-фотометрия по 6-му расширению `EPSF` для плотных полей (Crowded Fields)

- **Проблема:** Сейчас апертурная фотометрия в `app/services/spectrum_extractor.py` (`extract_sample`) суммирует поток в фиксированном круглом окне и оценивает фон по кольцу (annulus). Вблизи Млечного Пути или рядом с яркой звездой (например, `eps_indi_b` рядом с ярким компонентом Epsilon Indi A) в апертуру или фон попадают крылья соседней звезды, что искажает спектр и S/N слабого коричневого карлика.
- **Цель:** Использовать 6-е расширение L2 FITS-файлов SPHEREx (`EPSF` — эффективная PSF детектора) для PSF-подгонки и вычитания блендированных соседей.
- **План реализации:**
  1. Проверить/расширить `app/services/fits_stream.py` (или `spectrum_extractor.py`) для опционального чтения HDU `EPSF` (с корректным fallback на гауссову PSF или текущую апертурную фотометрию, если HDU `EPSF` отсутствует).
  2. Добавить в `app/services/spectrum_extractor.py` функцию PSF-фотометрии (взвешенный метод наименьших квадратов с маской флагов и весами $1/\sigma^2$ из `VARIANCE`), а при наличии рядом источников из `gaia_sources` (или `overlay`) в радиусе `< 12"` — совместную дебленд-подгонку нескольких PSF.
  3. Добавить в `SpectrumSummary` (`app/schemas/spectrum.py`) поле `photometry_method: Literal["aperture", "psf", "psf_deblended"] = "aperture"` и опциональное поле `crowding_warning: bool = False`.
  4. Написать unit-тесты в `tests/test_pipeline.py` с синтетическим плотным полем (две звезды на расстоянии ~2–3 пикселя), подтверждающие, что PSF-деблендинг восстанавливает спектр слабого компонента точнее обычной апертуры.
- **Промпт для агента:**
  > Прочитай `docs/AGENT_HANDOVER.md` (Задача 6), `app/services/fits_stream.py` и `app/services/spectrum_extractor.py`. Создай новую ветку от текущей и реализуй PSF-фотометрию с поддержкой 6-го расширения `EPSF` и деблендингом близких звёзд из Gaia-каталога, добавив `photometry_method` в `SpectrumSummary`. Покрой тестами на синтетическом плотном поле, закоммить и обнови `AGENT_HANDOVER.md` + `FRONTEND_GUIDE.md`.

---

### Задача 7 (Высокий приоритет): Эндпоинт кривых блеска `GET /api/v1/target/lightcurve`

- **Проблема:** В пресетах `PRESET_TARGETS` присутствуют переменные звёзды (`rw_aurigae` — звезда типа T Тельца, `r_draconis` — долгопериодическая переменная типа Миры), а в `DiffAnomaly` фиксируется лишь разница между двумя выбранными эпохами (`variable_brightening` / `variable_dimmed`). Астрофизикам и UI нужен полный временной ряд яркости (кривая блеска) по всем эпохам и диапазонам волн.
- **Цель:** Добавить эндпоинт `GET /api/v1/target/lightcurve`, возвращающий калиброванный поток (в мкЯн / uJy) в центре цели и опционально контрольной звезды сравнения в зависимости от времени `mjd` и детектора/полосы (`D1..D6` или группировка по длинам волн).
- **План реализации:**
  1. Добавить Pydantic-схемы `LightcurvePoint`, `LightcurveBandStats`, `LightcurveResponse` в `app/schemas/target.py` (или `app/schemas/lightcurve.py`).
  2. Реализовать метод `PipelineService.get_lightcurve(ra, dec, fov, band)` в `app/services/pipeline.py`, используя уже загруженные/закешированные FITS-кадры и апертурную/PSF-фотометрию по каждому отдельному кадру или мини-эпохе с указанием `mjd`, `detector`, `wavelength_um`, `flux_ujy`, `flux_err_ujy`, `snr`, `flagged`.
  3. Вычислить метрики переменности: `amplitude_ujy`, `variability_index` (например, $\chi^2_{\text{red}}$ относительно константы), `is_variable: bool`.
  4. Экспонировать роут `GET /api/v1/target/lightcurve` в `app/api/v1/targets.py` и добавить тест в `tests/test_pipeline.py`.
- **Промпт для агента:**
  > Прочитай `docs/AGENT_HANDOVER.md` (Задача 7), `app/api/v1/targets.py` и `app/services/pipeline.py`. В отдельной ветке реализуй эндпоинт `GET /api/v1/target/lightcurve` для построения кривых блеска по всем экспозициям SPHEREx (с разделением по детекторам `D1..D6` и метрикой переменности), добавь схемы, тест в `tests/test_pipeline.py`, документацию в `FRONTEND_GUIDE.md` и обнови `AGENT_HANDOVER.md`.

---

### Задача 8 (Средний приоритет): Пакетный сканер слепого поиска `POST /api/v1/survey/scan`

- **Проблема:** Сейчас пользователь должен заранее знать точные координаты `(ra, dec)` интересного объекта или выбирать один из 5 пресетов. Чтобы действительно «искать новое» на небе, нужен режим обзора области неба по сетке или списку координат с автоматическим ранжированием находок.
- **Цель:** Реализовать эндпоинт слепого поиска / пакетного сканирования поля и выдачу отсортированного рейтинга кандидатов (коричневых карликов и кандидатов в Планету X).
- **План реализации:**
  1. Спроектировать схемы `SurveyScanRequest` (центр области `ra, dec`, `radius_deg` или явный список координат, ограничения на число точек чтобы не перегрузить сервер) и `SurveyScanCandidate` / `SurveyScanResponse` с единым скорингом `candidate_score` (учитывающим отсутствие оптического двойника Gaia, высокий `spherex_snr`, красный цвет `ir_color_index` / WISE `W1 - W2`, `ch4_absorption_depth`, наличие `moving_dipole`, отсутствие совпадения со `known_minor_planets`, а также `motion_model == "parallax_dominated"`).
  2. Реализовать функцию расчёта интегрального `candidate_score` и оркестрацию сканера с ограничением конкурентности (`asyncio.Semaphore`), прогрессом и фоновым выполнением / кешированием результатов.
  3. Добавить роут `POST /api/v1/survey/scan` (и при необходимости `GET /api/v1/survey/scan/{job_id}`) и покрыть юнит-тестами скоринг + эндпоинт на мок-данных.
- **Промпт для агента:**
  > Прочитай `docs/AGENT_HANDOVER.md` (Задача 8) и `app/services/pipeline.py`. Создай ветку и реализуй сканер слепого поиска `POST /api/v1/survey/scan` с многофакторным скорингом кандидатов в Y-карлики / Планету X (учитывая Gaia, WISE, SkyBoT, CH4/Teff и `motion_model`), тестами в `tests/test_pipeline.py` и обновлением документации.

---

### Задача 9 (Средний приоритет): Экспорт научного архива `GET /api/v1/target/export`

- **Проблема:** Для публикации отчётов, работы в TOPCAT / SAOImageDS9 / Aladin или сохранения результатов хакатона пользователю и исследователю нужна кнопка скачивания сырых и обработанных данных в стандартных астрономических форматах.
- **Цель:** Добавить эндпоинт экспорта данных по цели в форматах `csv`, `votable` (IVOA XML), `fits` (многослойный FITS-куб / вырезки эпох) и `zip` (полный научный пакет: JSON-метаданные + спектр CSV/VOTable + PNG-слои).
- **План реализации:**
  1. Создать модуль `app/services/exporter.py` с генераторами CSV и IVOA VOTable (`astropy.io.votable` / `astropy.table.Table`), а также сборкой ZIP-архива (или FITS `HDUList` через `astropy.io.fits`) в памяти (`io.BytesIO` + `StreamingResponse`).
  2. Добавить эндпоинт `GET /api/v1/target/export` с query-параметрами `ra`, `dec`, `fov`, `band` и `format: Literal["csv", "votable", "fits", "zip"]` и корректными заголовками `Content-Disposition: attachment; filename=...`.
  3. Проверить на граничных случаях (цель с частичным кешем или отсутствующим спектром — возвращать понятный HTTP статус/описание) и добавить тесты валидности экспортируемых CSV/VOTable/FITS/ZIP в `tests/test_pipeline.py`.
- **Промпт для агента:**
  > Прочитай `docs/AGENT_HANDOVER.md` (Задача 9) и `app/api/v1/targets.py`. В отдельной ветке реализуй модуль `app/services/exporter.py` и эндпоинт `GET /api/v1/target/export` с поддержкой форматов `csv`, `votable`, `fits` и `zip`, добавь тесты валидации выгрузки и обнови `FRONTEND_GUIDE.md` и `AGENT_HANDOVER.md`.
