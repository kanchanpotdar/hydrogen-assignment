"""Compare Phase 3 presentation changes with the preserved V3 numeric fields."""
import csv
import io
from tools.reproduce_phase1 import compare


def compare_csv(filename, original, current):
    reader=csv.DictReader(io.StringIO(current.decode()))
    rows=list(reader)
    columns=list(reader.fieldnames)
    aliases={'signed_error':'abs_err', 'signed_difference':'abs_diff',
             'mean_species_RMSE':'species_curve_RMSE'}
    for row in rows:
        if 'absolute_error' in row:
            if float(row['absolute_error']) != abs(float(row['signed_error'])):
                raise AssertionError('absolute_error differs from abs(signed_error)')
        if filename=='smr_summary.csv':
            row['section']=row['section'].replace('(displayed Figure 4 segment)','(validation point)').replace(
                'P-R2 full 40 mm (experimental sampling outlet)',
                'P-R2 full 40 mm (extrapolation beyond paper Fig. 4)')
    stream=io.StringIO(newline='')
    selected=[c for c in columns if c!='absolute_error']
    writer=csv.DictWriter(stream,fieldnames=[aliases.get(c,c) for c in selected]);writer.writeheader()
    for row in rows:writer.writerow({aliases.get(c,c):row[c] for c in selected})
    result=compare(filename,original,stream.getvalue().encode())
    result['byte_equal']=original==current
    result['numeric_fields_exact']=all(v==0 for v in result['max_absolute_difference'].values())
    result['passed']=result['passed'] and result['numeric_fields_exact']
    result['allowed_presentation_mapping']=aliases
    return result
