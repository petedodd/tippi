# tippi
tbc

```
.
├── dataprep
│   ├── graphs
│   └── outdata
├── inference
│   ├── data
│   ├── graphs
│   ├── outdata
│   └── stan
└── model
    ├── data
    ├── graphs
    ├── indata
    └── outdata
```



# TODO

- identify bl data in report
- trace ASM & BC inputs
- finish commenting and tidying modeldata.R
- particular attention to including extra countries
- other York thresholds
- need to split out HIV entry ATT when cost merging?
- different screening costs?
- query zero screening costs
- CMR LSO missing cascade tabs
- NAs in results for ATT - down to missing countries in ASM <- trace
- understanding (HH vs facility CT) (also HIV vs not split)
- update WHO notifications
- consolidate other TODO s

test

Questions:
- ACF using FB screen cost?


# TODO from Nyasha call

- need to check unit costs, traceperhhcpt, totindexscreen - all mean the right thing
-> unit cost is per child contact screened
-> use data in Martina email to update HHCM cascade/split
- check logic of nhh ACF bits
- go back include the right entry-point splits

ACF:
- make sure that case-finding here does not include HHCM ACF
 -> query about whether INT row 2 includes HHCM screenings or not?
 (probably does include)
 (expect smaller; should be able to manage if not - see pt 3)
- HIV vs non-HIV entrypoint
-> what data has this in? new spreadsheet

Any unused costs? X-ray?
probably enters in analogous way to Xpert

see also joint intervention


Are FB vs CB splits same SoC vs INT?
---- 

PETE
- spreadsheet on splits for PT and ICF
- think about where Xray goes
- NA bug in ATT
- HH screen split (from ICF)
- overall cost drivers correct
- structuring the appendix
- trace ASM & BC inputs

NYASHA
- costing code base
- add to appendix on costing (relies on Pete outlining)
- York thresholds?
- fu w/Sushant


# New tables

## table 1

Simple logic:

| intervention | activity | country |
| ATT SOC      |          |         |
| ATT INT      |          |         |
| PT           |          |         |


Detailed items:

| intervention | activity                                   |
| ATT SOC      | Screened                                   |
| ATT SOC      | Presumptive TB                             |
| ATT SOC      | Tested with Xpert                          |
| ATT SOC      | TB diagnosed                               |
| ATT SOC      | TB treated                                 |
| ATT SOC      | Cost, $ (SD)                               |
| ATT INT      | Screened                                   |
| ATT INT      | Presumptive TB                             |
| ATT INT      | Tested with Xpert                          |
| ATT INT      | TB diagnosed                               |
| ATT INT      | TB treated                                 |
| ATT INT      | Cost, $ (SD)                               |
| TPT INT      | Ratio of PT initiations to ATT initiations |
| TPT INT      | PT via HIV clinic                          |
| TPT INT      | HHCM community based                       |
| TPT INT      | Started on PT                              |
| TPT INT      | Households screened                        |
| TPT INT      | Children screened                          |
| TPT INT      | Presumptive TB                             |
| TPT INT      | Diagnosed TB                               |
| TPT INT      | Cost, $ (SD)                               |

## table 2

Simple logic:

| internvetion | SOC | INT | Difference |
|  ICF         |     |     |            |
|  HHCM        |     |     |            |
|  combined    |     |     |            |


Detiled items

| Country | Started treatment | Cost | Started treatment | Cost | Started PT | Incident TB | Started ATT | TB deaths | Life-years (LYs) | Discounted LYs | Cost | ICER |
|         |                   |      |                   |      |            |             |             |           |                  |                |      |      |


To think: per TPT, per ATT, per...

Xray costs
