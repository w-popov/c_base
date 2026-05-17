#!/usr/bin/env python3

"""
ДЗ. 2 занятие. Системы счисления
Задачи №1 и №2.
Перевод чисел из десятичной системы в шестнадцатеричную и обратно.
"""
from enum import StrEnum

class TermColors(StrEnum):
    YE = '\033[1;33m'
    GR = '\033[1;32m'
    BL = '\033[1;36m'
    NC = '\033[0m'
    CLEAR_SCREEN = '\033[2J'


def dec_to_hex(dec: str) -> str:
    """Перевод десятичного числа в шестнадцатеричное методом деления на 16."""
   
    print(f"\n{TermColors.YE}Перевод {dec} в шестнадцатеричную систему:{TermColors.NC}")
    try:
        dec = int(dec)
    except ValueError:
        raise ValueError(f"Некорректное десятичное число: {dec}")

    remaiders = []
    str_solution = ""
    if not dec:
        return 0
    while dec:
        remaiders = [dec % 16] + remaiders
        str_solution += f"{dec} / 16 = {dec // 16} | Остаток: {dec % 16}\n"
        dec //= 16
    reply = "".join(format(x, "X") for x in remaiders)
    print(str_solution, f"{TermColors.BL}--> Ответ: {reply}{TermColors.NC}\n")



def hex_to_dec(hex_str: str):
    """Перевод шестнадцатеричного числа в десятичное схема Горнера."""
    try:
        int(hex_str, 16)
    except ValueError:
        raise ValueError(f"Некорректное шестнадцатеричное число: {hex_str}")

    print(f"{TermColors.GR}Перевод {hex_str} в десятичную систему:{TermColors.NC}")
    result = 0
    str_solution = ""
    for i, hex_digit in enumerate(hex_str):
        digit_value = int(hex_digit, 16)
        if not i:
            result = digit_value
        else:
            prev_res_ = result
            result = result * 16 + digit_value
            str_solution += f"{prev_res_} * 16 + {digit_value} = {result}\n"    
    print(str_solution, f"{TermColors.BL}--> Ответ: {result}{TermColors.NC}\n")
    


def main() -> None:
    print(TermColors.CLEAR_SCREEN)

    numbers = "12345678", "1000000",
    functions = dec_to_hex, hex_to_dec,
    
    for funcs in functions:
        for args in numbers:
            funcs(args)
    

if __name__ == "__main__":
    main()