#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int square_power (int number);                                      // B2
int summ_squares_range (int a, int b);                              // B3
const char *is_it_number_have_3_digit (unsigned number);            // B4
unsigned summ_digits_of_number (unsigned number);                   // B5
const char *is_both_equal_digits_in_number (int number);            // B6
const char *is_equal_digits_in_number (int number);                 // B7
const char *is_nine_digit_in_number (int number);                   // B8
const char *is_even_all_digits_in_number (int number);              // B9
const char *is_in_ascending_order (int number);                     // B10
unsigned reverse_number (unsigned number);                          // B11
int calculate_min_max_digits_number (unsigned *min, unsigned *max,
                                     unsigned _number);             // B12

int calculate_even_odd_digits_number (unsigned *even, unsigned *odd,
                                      unsigned _number);            // B13

unsigned number_of_even_numbers (char *, const int);                // B15

/* - - - - - - - B17 - - - - - - - - - - - - - - - - - - - -*/
struct UniquePtr_u
{
    unsigned *u_ptr;
};

void auto_free (struct UniquePtr_u *);

// Автоматическое освобождение памяти
#define UNIQE_PTR_U __attribute__((cleanup(auto_free))) struct UniquePtr_u

struct UniquePtr_u happy_numbers (unsigned);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const char *to_lowercase (char *text, char *lower,
                          const int size);                          // B21

#ifdef __cplusplus
}
#endif