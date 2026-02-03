# SARI data from Belo Horizonte

A dataset containing reports from Severe Acute Respiratory Illness
(SARI) from 2020 to April 2022 by week.

## Usage

``` r
noticeSARI
```

## Format

A data frame with 65404 rows and 7 variables:

- ref.week:

  The reference week, counting since the monitoring begun.

- reported.1.week:

  The number of cases occurred in the period and reported until the 1
  week after the reference week.

- reported.2.week:

  The number of cases occurred in the period and reported until the 2
  weeks after the reference week.

- reported.4.week:

  The number of cases occurred in the period and reported until the 4
  weeks after the reference week.

- reported.6.week:

  The number of cases occurred in the period and reported until the 6
  weeks after the reference week.

- reported.8.week:

  The number of cases occurred in the period and reported until the 8
  weeks after the reference week.

- reported.12.week:

  The number of cases occurred in the period and reported until the 12
  weeks after the reference week.

- occured:

  The total number of cases reported (at any time).

## Source

<https://datasus.saude.gov.br/informacoes-de-saude-tabnet/>
