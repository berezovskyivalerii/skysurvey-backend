# Документация API для фронтенд-команды (SkySurvey SPHEREx)

Интерактивный Swagger со всеми схемами и тестовыми кнопками: **`GET /docs`**  
Базовый префикс всех эндпоинтов: **`/api/v1`**  
Авторизация: **не требуется** (CORS открыт для всех доменов `*`).

---

## 1. Из каких блоков состоит интерфейс и откуда брать данные

Весь главный экран приложения строится вокруг **одного основного запроса** `GET /api/v1/target`, а остальные эндпоинты отвечают за интерактивные действия пользователя (ползунок яркости, клик по звезде для просмотра спектра, сравнение произвольных дат).

```text
+-----------------------------------------------------------------------------------+
| [Пресеты: WISE 0855 | Звезда Барнарда | ...]   RA: [133.7951] Dec: [-7.2451] [Go] |
+-----------------------------------------+-----------------------------------------+
|                                         |  ВЕРДИКТ (crossmatch.classification)    |
|         ГЛАВНЫЙ ВЬЮПОРТ 512x512         |  [ CANDIDATE / KNOWN STAR ]             |
|                                         |  Собств. движение, расстояние, SNR      |
|   Слои: (o) Ч/Б (image)                 +-----------------------------------------+
|         ( ) Цветной ИК (rgb)            |  ГРАФИК СПЕКТРА (0.75 - 5.0 мкм)        |
|         ( ) Карта движения (diff)       |  Данные: spectrum_data.chart_data       |
|         [x] Показать маску (flags)      |  Полосы H2O / CH4 / CO:                 |
|         [x] Кружки звезд Gaia (overlay) |  spectrum_data.absorption_bands_ref...  |
|         [x] Рамки движения (anomalies)  |                                         |
|                                         |                                         |
|  [ |> Play Blink ]  Дата: [2025-05]     |                                         |
|  Ползунок яркости (min_cut / max_cut)   |                                         |
+-----------------------------------------+-----------------------------------------+
```

---

## 2. Все эндпоинты и сценарии использования

### 2.1. `GET /api/v1/presets` — Кнопки демо-объектов в шапке
Возвращает 5 готовых объектов (коричневые карлики, летящая звезда Барнарда, переменные звезды).

**Запрос:**
```http
GET /api/v1/presets
```

**Ответ:**
```json
{
  "presets": [
    {
      "id": "wise_0855",
      "name": "WISE J085510.83-071442.5",
      "category": "brown_dwarf_candidate",
      "ra": 133.795129,
      "dec": -7.245147,
      "fov": 0.1,
      "band": "red",
      "description": "Coldest known brown dwarf (~250 K, Y4, d=2.28 pc, PM=8.1\"/yr)..."
    }
  ]
}
```
**Что делать на фронте:** Отрендерить кнопки быстрого выбора. При клике подставлять `ra`, `dec`, `fov` в поиск и вызывать `GET /api/v1/target`.

---

### 2.2. `GET /api/v1/target` — Главный запрос участка неба

**Параметры Query:**
| Параметр | Тип | По умолчанию | Описание |
| :--- | :--- | :--- | :--- |
| `ra` | `float` | *обязательный* | Прямое восхождение в градусах (`0 .. 360`), например `133.7951` |
| `dec` | `float` | *обязательный* | Склонение в градусах (`-90 .. +90`), например `-7.2451` |
| `fov` | `float` | `0.1` | Размер кадра в градусах (`0.05` — крупный зум, `0.1` — стандарт, `0.2` — широкий обзор) |
| `band` | `string` | `"all"` | Фильтр детекторов: `"all"`, `"blue"`, `"red"`, `"D1"`..`"D6"` |
| `sync` | `bool` | `true` | `true` — ждать готовности ответа; `false` — мгновенный ответ со статусом `processing` для прогресс-бара |

#### Два режима вызова (на выбор фронтенда):
1. **Простой режим (`sync=true`, по умолчанию):**
   - Если координаты уже есть в базе Supabase — ответ придет **за ~50 мс**.
   - Если координаты запрашиваются впервые — бэкенд скачает снимки из NASA, обработает их и вернет готовый JSON через **~10–20 секунд**.
2. **Режим с прогресс-баром (`sync=false`):**
   - Запрос `GET /api/v1/target?ra=...&dec=...&fov=0.1&sync=false` отвечает **мгновенно (HTTP 202)** и возвращает:
     ```json
     {
       "status": "processing",
       "progress": {
         "stage": "Resolving Datalink & downloading 18 cutouts",
         "percent": 35
       },
       "poll_after_ms": 1500
     }
     ```
   - Пока `status === "processing"`, фронт показывает прогресс-бар (`progress.percent` и `progress.stage`) и повторяет запрос каждые `poll_after_ms` (1.5 сек), пока `status` не станет `"ready"` (или `"partial"` — картинки уже готовы, можно рисовать!).

#### Структура ответа `GET /api/v1/target` (когда `status = "ready"`):

```json
{
  "target_key": "133.7951_-7.2451_0.100_all",
  "status": "ready",
  "ra": 133.7951,
  "dec": -7.2451,
  "fov": 0.1,
  "band": "all",
  "cache": { "hit": true, "stale": false },
  "n_observations": 258,

  "display": {
    "width": 512,
    "height": 512,
    "pixel_scale_arcsec": 0.7087,
    "wcs": {
      "CRVAL1": 133.7951,
      "CRVAL2": -7.2451,
      "CRPIX1": 256.5,
      "CRPIX2": 256.5,
      "CDELT1": -0.0001968587,
      "CDELT2": -0.0001968587
    },
    "stretch": {
      "stretch_type": "asinh",
      "vmin": -0.02478,
      "vmax": 0.10342,
      "unit": "MJy/sr"
    },
    "data_range": {
      "p0_5": -0.01853,
      "p50": -0.00006,
      "p99_5": 0.67659,
      "min": -0.30131,
      "max": 3.15295
    }
  },

  "epochs": [
    {
      "index": 0,
      "label": "2025-05",
      "date_obs": "2025-05-06T02:34:27Z",
      "mjd_mid": 60801.10726,
      "n_frames": 6,
      "detectors": [1, 2, 3, 4, 5, 6],
      "layers": {
        "image": "https://.../epoch_00_2025-05-06_image.png?v=1.0.0",
        "rgb": "https://.../epoch_00_2025-05-06_rgb.png?v=1.0.0",
        "variance": "https://.../epoch_00_2025-05-06_variance.png?v=1.0.0",
        "flags": "https://.../epoch_00_2025-05-06_flags.png?v=1.0.0"
      },
      "known_minor_planets": [
        {
          "number": "12345",
          "name": "2018 AB12",
          "object_type": "MBA",
          "ra": 133.7912,
          "dec": -7.2411,
          "v_mag": 18.4,
          "dra_arcsec_hr": -24.1,
          "ddec_arcsec_hr": 5.3,
          "separation_arcsec": 19.8,
          "x": 274.2,
          "y": 241.5
        }
      ]
    },
    {
      "index": 1,
      "label": "2025-11",
      "date_obs": "2025-11-07T09:05:37Z",
      "mjd_mid": 60986.37891,
      "n_frames": 6,
      "detectors": [1, 2, 3, 4, 5, 6],
      "layers": {
        "image": "https://.../epoch_01_2025-11-07_image.png?v=1.0.0",
        "rgb": "https://.../epoch_01_2025-11-07_rgb.png?v=1.0.0",
        "variance": "https://.../epoch_01_2025-11-07_variance.png?v=1.0.0",
        "flags": "https://.../epoch_01_2025-11-07_flags.png?v=1.0.0"
      }
    }
  ],

  "crossmatch": {
    "classification": "candidate",
    "reason": "5.9-sigma SPHEREx detection with no Gaia DR3 optical counterpart within 10\"",
    "catalog": "Gaia DR3",
    "match": null,
    "wise_match": {
      "source_id": "J085510.83-071442.5",
      "catalog": "CatWISE2020",
      "ra": 133.795125,
      "dec": -7.245139,
      "ra_obs": 133.772891,
      "dec_obs": -7.243211,
      "separation_arcsec": 0.65,
      "w1mpro": 17.82,
      "w2mpro": 13.89,
      "w1_w2_color": 3.93,
      "w1snr": 4.2,
      "w2snr": 28.5,
      "pmra_mas_yr": -8051.0,
      "pmdec_mas_yr": 698.0,
      "baseline_pm_arcsec_yr": 8.081,
      "is_ultracool_color": true
    },
    "spherex_snr": 5.89,
    "high_proper_motion": false,
    "overlay": [
      {
        "source_id": "5756955785394371328",
        "ra": 133.7970073,
        "dec": -7.25791,
        "g_mag": 13.29,
        "pm_total_mas_yr": 17.56,
        "x": 246.39,
        "y": 321.07
      }
    ]
  },

  "difference": {
    "epoch_a": 0,
    "epoch_b": 2,
    "baseline_days": 361.57,
    "diff_image_url": "https://.../diff_e00_e02.png?v=1.0.0",
    "deep_stack_url": "https://.../deep_stack.png?v=1.0.0",
    "max_dipole_snr": 8.4,
    "anomalies": [
      {
        "anomaly_type": "moving_dipole",
        "ra": 133.7949,
        "dec": -7.2452,
        "x": 257.4,
        "y": 256.8,
        "snr": 8.4,
        "x_from": 255.1,
        "y_from": 255.9,
        "shift_arcsec": 4.2,
        "pm_arcsec_yr": 4.24,
        "delta_flux_ujy": 145.6,
        "motion_model": "linear_proper_motion",
        "parallax_arcsec": null
      }
    ]
  },

  "spectrum_data": { "..." : "см. раздел 2.4" }
}
```

---

### Как фронтенду работать с блоками из `GET /api/v1/target`:

#### А. Таймлайн дат и анимация мигания (`epochs`)
- Число эпох (`epochs.length`) зависит от точки неба (обычно **2 или 3 эпохи**, например `["2025-08", "2026-03"]` или `["2025-05", "2025-11", "2026-05"]`). Не хардкодьте число 3 — рендерите кнопки дат динамически по массиву `epochs`.
- **Анимация «Blink» (как в WiseView):** запустите `setInterval` на `400–600 мс`, который переключает текущий индекс эпохи `activeEpoch = (activeEpoch + 1) % epochs.length`.
- Все картинки уже выровнены по координатам и яркости. При переключении `epochs[activeEpoch].layers.image` (или `layers.rgb`) неподвижные звезды стоят на месте, а движущийся объект прыгает!
- **Слой `layers.flags`:** это полупрозрачный PNG (`512×512`) с оранжевыми/красными пикселями там, где попали космические лучи или битые пиксели матрицы. Его можно класть **поверх** основной картинки (`position: absolute; top: 0; left: 0; pointer-events: none`) по чекбоксу «Показать артефакты/битые пиксели».

#### Б. Бейдж классификации и наложение известных звезд (`crossmatch`) и астероидов (`known_minor_planets`)
- Поле `crossmatch.classification` принимает значения:
  - **`"candidate"`** — **ГЛАВНАЯ НАХОДКА!** Объект ярко светится в инфракрасном телескопе SPHEREx, но отсутствует в оптическом каталоге Gaia DR3. Подсветите ярким бейджем **«Кандидат (Коричневый карлик / Планета X)»**.
  - **`"known_star"`** — объект совпал с известной звездой Gaia DR3 (в `crossmatch.match` лежат её ID, расстояние в парсеках `distance_pc`, блеск `phot_g_mean_mag` и скорость движения `pm_total_mas_yr`). Если `crossmatch.high_proper_motion === true`, покажите бейдж **«Быстролетящая звезда (High Proper Motion)»**.
  - **`"no_detection"`** — в центре кадра только фоновый шум.
- **`crossmatch.wise_match`** (если найден ИК-источник в `CatWISE2020` или `AllWISE`): содержит `w1mpro` (3.35 мкм), `w2mpro` (4.60 мкм), цвет `w1_w2_color` (`W1 - W2`), флаг `is_ultracool_color` (`true` при `W1 - W2 >= 1.5 mag`, типично для ультрахолодных T/Y-карликов), а также `baseline_pm_arcsec_yr` — собственное движение за ~10 лет между эпохой WISE и SPHEREx (и перенесенные координаты `ra_obs`, `dec_obs`). Выведите отдельную карточку **«WISE 10-yr Cross-Check»**.
- Массив `crossmatch.overlay`: содержит уже готовые пиксельные координаты **`(x, y)` в диапазоне `0..512`**. Нарисуйте в `<svg viewBox="0 0 512 512">` поверх картинки тонкие голубые кружки в точках `(item.x, item.y)` — пользователь сразу увидит, какие источники на снимке уже известны науке, а какие — нет!
- Массив `epochs[i].known_minor_planets`: список известных малых тел Солнечной системы (астероиды/кометы из IMCCE SkyBoT на момент `mjd_mid` каждой эпохи) с готовыми координатами `(x, y)` в `0..512`, `name`, `object_type`, `v_mag`, `dra_arcsec_hr`, `ddec_arcsec_hr`. Нарисуйте их оранжевыми ромбами или пунктирными кружками поверх кадра текущей эпохи, чтобы не перепутать транзитный астероид с Планетой X!

#### В. Карта движения и авто-поиск аномалий (`difference`)
- `difference.diff_image_url` — картинка разности между самой первой и самой последней датой (`512×512 PNG`). Неподвижные звезды на ней вычтены в белый/нейтральный фон; синее пятно — откуда объект ушел, красное — куда пришел.
- `difference.anomalies` — массив автоматически найденных движущихся объектов (`"moving_dipole"`) и вспыхнувших/потускневших переменных звезд (`"variable_brightened"` / `"variable_dimmed"`):
  - Для `"moving_dipole"` нарисуйте поверх картинки стрелочку или рамку от `(x_from, y_from)` к `(x, y)` и подпишите скорость `pm_arcsec_yr` (`"/год`).
  - Поля траекторной модели (при наличии 3+ эпох): `motion_model` (`"linear_proper_motion"` — линейное движение близкого коричневого карлика / быстрой звезды; `"parallax_dominated"` — преобладает годовой параллакс Земли, характерный для далёкого тела Солнечной системы вроде Планеты X; `"ambiguous"` — недостаточно эпох или неоднозначно) и `parallax_arcsec` (оценка амплитуды параллакса в угловых секундах). Покажите соответствующий бейдж рядом с диполем!

---

### 2.3. `POST /api/v1/target/restretch` — Ползунки яркости и контраста

Когда пользователь хочет «вытянуть» тусклый объект или убрать засветку фона, он крутит ползунки нижней (`min_cut`) и верхней (`max_cut`) границы яркости или меняет алгоритм (`asinh`, `log`, `zscale`).

**Откуда взять стартовые значения ползунков (`min`, `max`, `step`):**
Из ответа `GET /api/v1/target` в блоке `display`:
- Текущие значения ползунков по умолчанию:
  - `min_cut = display.stretch.vmin` (например, `-0.024`)
  - `max_cut = display.stretch.vmax` (например, `0.103`)
- Диапазон, в котором ходят ползунки (min/max самого `<input type="range">`):
  - Нижний предел слайдера: `display.data_range.p0_5` (например, `-0.05`)
  - Верхний предел слайдера: `display.data_range.p99_5` (например, `0.70`)

**Запрос (`POST /api/v1/target/restretch`):**
```json
{
  "ra": 269.4485,
  "dec": 4.7394,
  "fov": 0.1,
  "band": "all",
  "stretch_type": "asinh",
  "min_cut": -0.02,
  "max_cut": 0.25,
  "asinh_a": 0.08
}
```
*(Важно: `max_cut` всегда должен быть строго больше `min_cut`. Если передать `"min_cut": null, "max_cut": null`, бэкенд сам подберет авто-границы).*

**Ответ (~15–30 мс):**
```json
{
  "target_key": "269.4485_+4.7394_0.100_all",
  "stretch": {
    "stretch_type": "asinh",
    "vmin": -0.02,
    "vmax": 0.25,
    "unit": "MJy/sr"
  },
  "width": 512,
  "height": 512,
  "epochs": [
    { "index": 0, "image": "data:image/png;base64,iVBORw0KGgoAAA..." },
    { "index": 1, "image": "data:image/png;base64,iVBORw0KGgoAAA..." }
  ],
  "source": "local_cache",
  "elapsed_ms": 14.8
}
```
**Что делать на фронте:**
- При движении ползунка яркости (с `debounce` `100–150 мс`) отправлять `POST /api/v1/target/restretch`.
- В ответе приходят обновленные картинки **сразу для всех дат (`epochs`)** в виде `data:image/png;base64,...`.
- Просто замените `layers.image` у каждой эпохи на полученный `epochs[i].image` — тогда и текущий кадр, и анимация мигания (**Blink**) между датами сразу продолжат работать с новым контрастом без перезагрузки файлов по сети!

---

### 2.4. `GET /api/v1/target/spectrum` — График спектра (102 канала SPHEREx)

Спектр центрального объекта **уже лежит** внутри главного ответа `GET /api/v1/target` в поле `spectrum_data`.
Отдельный эндпоинт `GET /api/v1/target/spectrum?ra=...&dec=...` нужен тогда, когда:
1. Пользователь **кликнул мышкой в любую точку или звезду на снимке**, и нужно быстро построить спектр именно для той точки, куда он кликнул.
2. Или нужно обновить только график спектра, не запрашивая картинки заново.

**Запрос:**
```http
GET /api/v1/target/spectrum?ra=133.7951&dec=-7.2451
```

**Ответ:**
```json
{
  "target": { "ra": 133.7951, "dec": -7.2451 },
  "spectrum_summary": {
    "total_channels": 102,
    "valid_points": 18,
    "flagged_points": 0,
    "n_exposures_used": 18,
    "mean_snr": 1.65,
    "max_snr": 3.94,
    "peak_wavelength_um": 4.2882,
    "peak_flux_ujy": 202.793,
    "ir_color_index": 1.95,
    "spectral_hint": "red_stellar_or_dusty_source",
    "flux_unit": "uJy",
    "estimated_temp_k": 310.0,
    "ch4_absorption_depth": 0.42,
    "h2o_absorption_depth": 0.31,
    "model_curve": [
      {"wavelength_um": 0.75, "model_flux": 0.12},
      {"wavelength_um": 3.30, "model_flux": 48.5},
      {"wavelength_um": 4.60, "model_flux": 198.4}
    ]
  },
  "chart_data": [
    {
      "channel_id": 6,
      "wavelength_um": 0.8408,
      "bandwidth_um": 0.021,
      "band": "Band 1",
      "detector": 1,
      "flux": 93.363,
      "flux_error": 56.845,
      "snr": 1.64,
      "is_flagged": false,
      "n_samples": 1
    }
  ],
  "absorption_bands_reference": [
    { "molecule": "H2O", "range_um": [1.1, 1.2], "label": "H2O (J-band)" },
    { "molecule": "CH4", "range_um": [1.6, 1.7], "label": "CH4 (H-band T/Y dwarf)" },
    { "molecule": "CH4", "range_um": [3.2, 3.4], "label": "CH4 (3.3 um fundamental)" },
    { "molecule": "CO2", "range_um": [4.2, 4.35], "label": "CO2 (4.26 um)" },
    { "molecule": "CO",  "range_um": [4.55, 4.8], "label": "CO (4.67 um fundamental)" }
  ]
}
```

**Как рисовать график на фронте (`Recharts` / `Chart.js`):**
- **Ось X:** `wavelength_um` (длина волны в мкм, от `0.75` до `5.0`, массив уже строго отсортирован по возрастанию).
- **Ось Y:** `flux` (поток в микроянских `μJy`).
- **Усы погрешности (Error Bars):** `flux_error` ($\pm 1\sigma$).
- **Флаги качества (`is_flagged`):** точки с `is_flagged: true` можно рисовать серым/прозрачным кружком или скрывать по чекбоксу «Скрыть шумные точки».
- **Вертикальные полосы молекул (`absorption_bands_reference`):** нарисуйте полупрозрачные вертикальные зоны (`<ReferenceArea x1={range_um[0]} x2={range_um[1]} />`) для воды ($\text{H}_2\text{O}$), метана ($\text{CH}_4$) и $\text{CO}$ — провалы графика внутри этих полос наглядно доказывают наличие холодной атмосферы у коричневого карлика!
- **Карточки метрик над графиком (`spectrum_summary`):**
  - `ir_color_index`: отношение яркости на `3–5 мкм` к `0.75–1.5 мкм` (если `> 2.5` — объект аномально красный/холодный).
  - `estimated_temp_k`: автоматическая оценка эффективной температуры чёрного тела Планка ($T_{\text{eff}}$ в Кельвинах, диапазон `200..6000` K; может быть `null`, если точек недостаточно).
  - `ch4_absorption_depth`: глубина полосы поглощения метана $\text{CH}_4$ на `3.2–3.4 мкм` (`0.0..1.0`, где `> 0.25` указывает на холодную атмосферу T/Y-карлика).
  - `h2o_absorption_depth`: глубина полосы поглощения водяного пара $\text{H}_2\text{O}$ на `2.55–2.85 мкм` (`0.0..1.0`).
  - `model_curve`: массив точек `{wavelength_um, model_flux}` для отрисовки плавной теоретической кривой (пунктирная линия поверх измеренных точек `chart`).
  - `spectral_hint`: человекочитаемый класс спектра:
    - `"ultra_cool_y_dwarf_or_planet_x_candidate"` — ультрахолодный Y-карлик / кандидат в Планету X (яркий пик на `4–5 мкм`).
    - `"cold_brown_dwarf_candidate"` — холодный коричневый карлик (T/L-класс с метановым поглощением).
    - `"red_stellar_or_dusty_source"` — красный карлик или пылевой источник.
    - `"normal_stellar_source"` — обычная горячая звезда (спектр падает к длинным волнам).

---

### 2.5. `GET /api/v1/target/diff` — Сравнение любых двух выбранных дат

В главном ответе `GET /api/v1/target` поле `difference` всегда сравнивает **первую (`0`)** и **последнюю** эпоху. Если у объекта 3 эпохи (`0`, `1`, `2`) и пользователь хочет сравнить, например, `0` и `1` или `1` и `2`:

**Запрос:**
```http
GET /api/v1/target/diff?ra=133.7951&dec=-7.2451&fov=0.1&epoch_a=0&epoch_b=1
```
**Ответ:** Тот же объект `DifferenceResult` (`diff_image_url`, `baseline_days`, `anomalies`).

---

## 3. Готовые рецепты кода для фронтенда (TypeScript / React)

### Рецепт 1: Как по клику мыши `(x, y)` на картинке `512×512` узнать небесные координаты `(ra, dec)` и запросить спектр этой точки

Бэкенд отдает в `display.wcs` точную привязку пикселей `512×512` к небу. Перевод клика мыши в `(ra, dec)` занимает 6 строчек кода без всяких астрономических библиотек:

```ts
export function pixelToRaDec(
  pixelX: number, // 0..512 (слева направо)
  pixelY: number, // 0..512 (сверху вниз, как в DOM mouse event)
  wcs: Record<string, number | string>
): { ra: number; dec: number } {
  const crval1 = Number(wcs.CRVAL1);
  const crval2 = Number(wcs.CRVAL2);
  const crpix1 = Number(wcs.CRPIX1); // 256.5
  const crpix2 = Number(wcs.CRPIX2); // 256.5
  const cdelt1 = Number(wcs.CDELT1); // шаг по RA (отрицательный, т.к. Восток слева)
  const cdelt2 = Number(wcs.CDELT2); // шаг по Dec

  const dx = (pixelX + 1.0) - crpix1;
  const dy = (pixelY + 1.0) - crpix2;

  const dec = crval2 + dy * cdelt2;
  const cosDec = Math.cos((crval2 * Math.PI) / 180.0) || 1.0;
  const ra = (crval1 + (dx * cdelt1) / cosDec + 360.0) % 360.0;

  return {
    ra: Number(ra.toFixed(6)),
    dec: Number(dec.toFixed(6)),
  };
}
```
**Как использовать при клике по картинке:**
1. В `onClick` на контейнере картинки `512×512` берете `x` и `y` клика.
2. Вызываете `const { ra, dec } = pixelToRaDec(x, y, targetData.display.wcs)`.
3. Вызываете `fetch(\`/api/v1/target/spectrum?ra=\${ra}&dec=\${dec}\`)` и обновляете график спектра!

---

### Рецепт 2: Как сделать плавный ползунок контраста (`/restretch`) с `debounce`

```ts
import { useState, useEffect } from "react";

export function useRestretch(
  ra: number,
  dec: number,
  fov: number,
  stretchType: "asinh" | "log" | "zscale",
  minCut: number | null,
  maxCut: number | null,
  onImagesUpdated: (epochImages: Record<number, string>) => void
) {
  useEffect(() => {
    if (minCut !== null && maxCut !== null && maxCut <= minCut) return;

    const timer = setTimeout(async () => {
      const res = await fetch("/api/v1/target/restretch", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          ra,
          dec,
          fov,
          stretch_type: stretchType,
          min_cut: minCut,
          max_cut: maxCut,
        }),
      });
      if (!res.ok) return;
      const data = await res.json();
      // Превращаем массив [{ index: 0, image: "data:image/png;base64,..." }] в словарь по индексам эпох
      const map: Record<number, string> = {};
      for (const ep of data.epochs) {
        map[ep.index] = ep.image;
      }
      onImagesUpdated(map);
    }, 120); // 120 мс debounce

    return () => clearTimeout(timer);
  }, [ra, dec, fov, stretchType, minCut, maxCut]);
}
```

---

### Рецепт 3: Загрузка `GET /api/v1/target` с поддержкой прогресс-бара (`sync=false`)

```ts
export async function fetchTargetWithProgress(
  ra: number,
  dec: number,
  fov: number,
  onProgress?: (stage: string, percent: number) => void
) {
  while (true) {
    const res = await fetch(
      `/api/v1/target?ra=${ra}&dec=${dec}&fov=${fov}&sync=false`
    );
    if (!res.ok) {
      throw new Error(`API error: ${res.status}`);
    }
    const data = await res.json();

    if (data.progress && onProgress) {
      onProgress(data.progress.stage, data.progress.percent);
    }

    // Как только статус стал "ready" или "partial" — все картинки и спектр уже в ответе!
    if (data.status === "ready" || data.status === "partial" || data.status === "no_data") {
      return data;
    }

    await new Promise((r) => setTimeout(r, data.poll_after_ms ?? 1500));
  }
}
```
