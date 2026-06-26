# Author - Eva Sinha, Pacific Northwest National Lab
#
# era5_hdd_cdd_to_gcam_xml.py
#
# Converts ERA5-derived HDD and CDD CSVs (output of era5_hdd_cdd_gcam_regions.py)
# into GCAM XML input files that follow the same structure as:
#   HDDCDD_constdd_no_GCM.xml   (32 global GCAM regions)
#   HDDCDD_constdds_USA.xml     (51 US states)
#
# APPROACH:
#   The ERA5 CSVs contain one row per region x GLU (or state x GLU).
#   These are collapsed to one value per region (or state) using a
#   land-area-weighted mean across GLUs.  The result is rounded to the
#   nearest integer to match the existing XML format.
#
# ERA5 PERIOD LABEL → XML HISTORICAL YEAR MAPPING:
#   ERA5 1980 (avg 1976–1980)  → XML year 1975
#   ERA5 1990 (avg 1986–1990)  → XML year 1990
#   ERA5 2005 (avg 2001–2005)  → XML year 2005
#   ERA5 2010 (avg 2006–2010)  → XML year 2010
#   ERA5 2015 (avg 2011–2015)  → XML year 2015
#   Future years 2020–2100     → hold 2015 value constant
#
# WORKFLOW:
#   Step 1 - Load and aggregate global region CSVs
#     For each of the 5 ERA5 labels, read gcam_regions_HDD_CDD_<label>.csv
#     and compute land_area_km2-weighted mean HDD and CDD per gcam_reg_name.
#
#   Step 2 - Load and aggregate USA state CSVs
#     For each of the 5 ERA5 labels, read gcam_usa_states_HDD_CDD_<label>.csv,
#     filter to state rows, and compute weighted mean per state abbreviation.
#
#   Step 3 - Build and write global XML (HDDCDD_ERA5_no_GCM.xml)
#     One <region> per GCAM region; comm and resid consumers share the same
#     degree-day values.
#
#   Step 4 - Build and write USA states XML (HDDCDD_ERA5_USA.xml)
#     One <region> per US state abbreviation.
#
# INPUT FILES (./hdd_cdd_outfiles/):
#   gcam_regions_HDD_CDD_1980.csv, ..., gcam_regions_HDD_CDD_2015.csv
#   gcam_usa_states_HDD_CDD_1980.csv, ..., gcam_usa_states_HDD_CDD_2015.csv
#
# OUTPUT FILES (./hdd_cdd_outfiles/):
#   HDDCDD_ERA5_no_GCM.xml   -- global GCAM regions
#   HDDCDD_ERA5_USA.xml      -- US states

import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HDD_CDD_DIR = Path('./hdd_cdd_outfiles')
OUT_DIR     = HDD_CDD_DIR

# ERA5 period label → GCAM XML historical year
ERA5_TO_XML_YEAR = {
    1980: 1975,
    1990: 1990,
    2005: 2005,
    2010: 2010,
    2015: 2015,
}
XML_HIST_YEARS = sorted(ERA5_TO_XML_YEAR.values())   # [1975, 1990, 2005, 2010, 2015]
FUTURE_YEARS   = list(range(2020, 2105, 5))           # 2020, 2025, ..., 2100

# US state abbreviations used to filter the USA CSV (which also contains
# global-region rows).
US_STATE_ABBREVS = {
    'AK', 'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA', 'HI',
    'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME', 'MI', 'MN',
    'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM', 'NV', 'NY', 'OH',
    'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'VT', 'WA',
    'WI', 'WV', 'WY',
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_and_aggregate(csv_path: Path, filter_values: set | None = None) -> pd.DataFrame:
    """Return land-area-weighted mean HDD and CDD per gcam_reg_name (integer).

    Parameters
    ----------
    csv_path      : Path to a gcam_regions_HDD_CDD or gcam_usa_states_HDD_CDD CSV.
    filter_values : If given, keep only rows whose gcam_reg_name is in this set.

    Returns
    -------
    DataFrame indexed by gcam_reg_name with columns HDD and CDD (int).
    """
    df = pd.read_csv(csv_path)
    if filter_values is not None:
        df = df[df['gcam_reg_name'].isin(filter_values)]

    df = df[df['land_area_km2'] > 0].copy()
    df['hdd_w'] = df['HDD'] * df['land_area_km2']
    df['cdd_w'] = df['CDD'] * df['land_area_km2']

    agg = df.groupby('gcam_reg_name').agg(
        total_land=('land_area_km2', 'sum'),
        hdd_sum   =('hdd_w',        'sum'),
        cdd_sum   =('cdd_w',        'sum'),
    ).reset_index()

    agg['HDD'] = (agg['hdd_sum'] / agg['total_land']).round().astype(int)
    agg['CDD'] = (agg['cdd_sum'] / agg['total_land']).round().astype(int)
    return agg.set_index('gcam_reg_name')[['HDD', 'CDD']]


def add_region_to_world(world_el: ET.Element, region_name: str,
                        year_vals: dict[int, dict[str, int]]) -> None:
    """Append a fully populated <region> element to a <world> element.

    Parameters
    ----------
    world_el    : Parent <world> XML element.
    region_name : Value for the 'name' attribute of <region>.
    year_vals   : {xml_year: {'HDD': int, 'CDD': int}} for all historical years.
                  Future years automatically repeat the 2015 value.
    """
    region_el = ET.SubElement(world_el, 'region', name=region_name)

    for consumer_name in ['comm', 'resid']:
        consumer  = ET.SubElement(region_el, 'gcam-consumer', name=consumer_name)
        node_in   = ET.SubElement(consumer, 'nodeInput', name=consumer_name)
        bldg_in   = ET.SubElement(node_in, 'building-node-input',
                                   name=f'{consumer_name}_building')

        for service_suffix, key in [('cooling', 'CDD'), ('heating', 'HDD')]:
            thermal = ET.SubElement(bldg_in, 'thermal-building-service-input',
                                    name=f'{consumer_name} {service_suffix}')

            for xml_year in XML_HIST_YEARS:
                el      = ET.SubElement(thermal, 'degree-days', year=str(xml_year))
                el.text = str(year_vals[xml_year][key])

            val_2015 = year_vals[2015][key]
            for yr in FUTURE_YEARS:
                el      = ET.SubElement(thermal, 'degree-days', year=str(yr))
                el.text = str(val_2015)


def write_xml(root_el: ET.Element, out_path: Path) -> None:
    """Indent and write an XML element tree with UTF-8 declaration."""
    ET.indent(root_el, space='    ')
    ET.ElementTree(root_el).write(out_path, xml_declaration=True, encoding='UTF-8')
    print(f"  Written: {out_path}")


def build_xml(by_year: dict[int, pd.DataFrame], label: str) -> ET.Element:
    """Build a <scenario> XML tree from per-year aggregated DataFrames.

    Parameters
    ----------
    by_year : {xml_year: DataFrame indexed by region_name with HDD/CDD columns}
    label   : Human-readable label for progress messages (e.g. 'global', 'USA states')
    """
    # Collect all region names present in every year (inner intersection).
    all_regions = sorted(
        set.intersection(*[set(df.index) for df in by_year.values()])
    )
    missing_any = set.union(*[set(df.index) for df in by_year.values()]) - set(all_regions)
    if missing_any:
        print(f"  WARNING ({label}): {len(missing_any)} region(s) absent from at "
              f"least one period and excluded: {sorted(missing_any)}")

    scenario = ET.Element('scenario')
    world    = ET.SubElement(scenario, 'world')

    for region_name in all_regions:
        year_vals = {
            yr: {'HDD': int(by_year[yr].loc[region_name, 'HDD']),
                 'CDD': int(by_year[yr].loc[region_name, 'CDD'])}
            for yr in XML_HIST_YEARS
        }
        add_region_to_world(world, region_name, year_vals)

    print(f"  Built XML for {len(all_regions)} {label} regions.")
    return scenario

# ---------------------------------------------------------------------------
# Step 1 – Load and aggregate global region CSVs
# ---------------------------------------------------------------------------
print(f"Starting ERA5 HDD/CDD → GCAM XML conversion at: {datetime.now()}\n")
print("Step 1: Loading global region HDD/CDD CSVs ...")

global_by_year: dict[int, pd.DataFrame] = {}
for era5_label, xml_year in sorted(ERA5_TO_XML_YEAR.items()):
    csv_path = HDD_CDD_DIR / f'gcam_regions_HDD_CDD_{era5_label}.csv'
    global_by_year[xml_year] = load_and_aggregate(csv_path)
    print(f"  ERA5 {era5_label} → XML {xml_year}: "
          f"{len(global_by_year[xml_year])} regions  "
          f"(HDD {global_by_year[xml_year]['HDD'].min()}–{global_by_year[xml_year]['HDD'].max()}, "
          f"CDD {global_by_year[xml_year]['CDD'].min()}–{global_by_year[xml_year]['CDD'].max()})")

# ---------------------------------------------------------------------------
# Step 2 – Load and aggregate USA state CSVs
# ---------------------------------------------------------------------------
print("\nStep 2: Loading USA state HDD/CDD CSVs ...")

state_by_year: dict[int, pd.DataFrame] = {}
for era5_label, xml_year in sorted(ERA5_TO_XML_YEAR.items()):
    csv_path = HDD_CDD_DIR / f'gcam_usa_states_HDD_CDD_{era5_label}.csv'
    state_by_year[xml_year] = load_and_aggregate(csv_path, filter_values=US_STATE_ABBREVS)
    print(f"  ERA5 {era5_label} → XML {xml_year}: "
          f"{len(state_by_year[xml_year])} states  "
          f"(HDD {state_by_year[xml_year]['HDD'].min()}–{state_by_year[xml_year]['HDD'].max()}, "
          f"CDD {state_by_year[xml_year]['CDD'].min()}–{state_by_year[xml_year]['CDD'].max()})")

# ---------------------------------------------------------------------------
# Step 3 – Build and write global XML
# ---------------------------------------------------------------------------
print("\nStep 3: Building global regions XML ...")
scenario_global = build_xml(global_by_year, label='global')
write_xml(scenario_global, OUT_DIR / 'HDDCDD_ERA5_no_GCM.xml')

# ---------------------------------------------------------------------------
# Step 4 – Build and write USA states XML
# ---------------------------------------------------------------------------
print("\nStep 4: Building USA states XML ...")
scenario_usa = build_xml(state_by_year, label='USA states')
write_xml(scenario_usa, OUT_DIR / 'HDDCDD_ERA5_USA.xml')

print(f"\nFinished at: {datetime.now()}")
