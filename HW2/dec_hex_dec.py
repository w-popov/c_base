#!/usr/bin/env python3

"""
ДЗ. 2 занятие. Системы счисления
Задачи №1 и №2.
Перевод чисел из десятичной системы в шестнадцатеричную и обратно.
"""

# первое число
NUMBER_1 = "12345678"

# второе число
NUMBER_2 = "1000000"


def dec_to_hex(dec: str) -> str:
    """Перевод десятичного числа в шестнадцатеричное методом деления на 16."""
   
    print(f"\n\033[1;33mПеревод {dec} в шестнадцатеричную систему:\033[0m")
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
    print(str_solution)
    return "".join(format(x, "X") for x in remaiders)



def hex_to_dec(hex: str) -> int:
    """Перевод шестнадцатеричного числа в десятичное схема Горнера."""
    try:
        int(hex, 16)
    except ValueError:
        raise ValueError(f"Некорректное шестнадцатеричное число: {hex_str}")

    print(f"\033[1;32mПеревод {hex} в десятичную систему:\033[0m")
    result = 0
    str_solution = ""
    for i, hex_digit in enumerate(hex):
        digit_value = int(hex_digit, 16)
        if not i:
            result = digit_value
        else:
            prev_res_ = result
            result = result * 16 + digit_value
            str_solution += f"{prev_res_} * 16 + {digit_value} = {result}\n"    
    print(str_solution)
    return result


def main() -> None:
    print("\033[2J") # очистить весь экран

    print(f"\033[1;36m--> Ответ: {dec_to_hex(NUMBER_1)}\033[0m\n")
    print(f"\033[1;36m--> Ответ: {dec_to_hex(NUMBER_2)}\033[0m\n")
    print("- " * 25, "\n")
    print(f"\033[1;36m--> Ответ: {hex_to_dec(NUMBER_1)}\033[0m\n")
    print(f"\033[1;36m--> Ответ: {hex_to_dec(NUMBER_2)}\033[0m\n")
    

if __name__ == "__main__":
    main()