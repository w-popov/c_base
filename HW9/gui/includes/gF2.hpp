#pragma once

#include <FL/Fl_Double_Window.H>
#include <vector>

class Fl_Multiline_Output;
class Fl_Input;
class Fl_Output;
class Fl_Button;
class Fl_Grid;
class Fl_Box;

constexpr const char *INFO2 =
    "ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г\n\
F2: Написать только одну функцию, которая ставит в начало массива все четные элементы, а в конец – все нечетные \
void sort_even_odd(int n, int a[]). \
Не нарушайте порядок следования чисел между собой.\n";

class WindowF2 : public Fl_Double_Window
{
  private:
    Fl_Input *in{nullptr};
    Fl_Output *out{nullptr};
    Fl_Button *btn{nullptr};
    Fl_Grid *grid{nullptr};
    Fl_Box *ibox{nullptr};
    Fl_Box *title{nullptr};
    Fl_Box *info_box{nullptr};
    std::vector<int> array;

  private:
    static void button_callback (Fl_Widget *, void *data);
    std::vector<int> &parse_input (const char *input);
    void print_output (const std::vector<int> &array);

  public:
    WindowF2(int w, int h, const char *title = "ДЗ-9: F2");
    void sort_eodd (int *array, int size);
};