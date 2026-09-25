
# 分组函数（从 Cell 10 复制）
def assign_group_A(row):
    c = row['Cyt_any']
    m = row['Mem_any']
    n = row['Nuc_any']
    if c == 1 and m == 0 and n == 0:
        return "C-only"
    elif c == 0 and m == 1 and n == 0:
        return "M-only"
    elif c == 0 and m == 0 and n == 1:
        return "N-only"
    elif c == 1 and m == 1 and n == 0:
        return "C&M"
    elif c == 1 and m == 0 and n == 1:
        return "C&N"
    elif c == 0 and m == 1 and n == 1:
        return "M&N"
    elif c == 1 and m == 1 and n == 1:
        return "C&M&N"
    else:
        return "None"

def get_state_at_time(row, time_label):
    c = row[f'{time_label}_Cyt']
    m = row[f'{time_label}_Mem']
    n = row[f'{time_label}_Nuc']
    if c == 1 and m == 0 and n == 0:
        return "C"
    elif c == 0 and m == 1 and n == 0:
        return "M"
    elif c == 0 and m == 0 and n == 1:
        return "N"
    elif c == 1 and m == 1 and n == 0:
        return "CM"
    elif c == 1 and m == 0 and n == 1:
        return "CN"
    elif c == 0 and m == 1 and n == 1:
        return "MN"
    elif c == 1 and m == 1 and n == 1:
        return "CMN"
    else:
        return "None"
