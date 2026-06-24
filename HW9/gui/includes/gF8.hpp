#pragma once

#include <FL/Fl_Double_Window.H>
#include <vector>
#include <string>

#ifdef _WIN32
#include <windows.h>
#include <FL/x.H>
void seticon(Fl_Window *win);
#endif

class F8
{
  private:
    std::vector<int> array;
  public:
    ~F8() = default;
    F8() = default;
    F8(const char *arr_values);
    F8(const F8&);
    F8(F8&&) noexcept;
    F8 &operator = (const F8& f8)
    {
        if (this != &f8)
        {
            array = f8.array;
        }
        return *this;
    }

    F8 &operator = (F8&& f8) noexcept
    {
        if (this != &f8)
        {
            array = std::move(f8.array);
        }
        return *this;
    }

    int missing_num();
    std::vector<int> &init(const char *input);
    bool is_empty() const;
    const std::vector<int>& get_array() const;

  private:
    std::vector<int> parse_input(const char *input);
};

class Fl_Scroll;
class Fl_Input;
class Fl_Output;
class Fl_Button;
class Fl_Grid;
class Fl_Box;

class WindowF8 : public Fl_Double_Window
{
  private:
    Fl_Grid *cells_grid{nullptr};
    Fl_Scroll *scroll{nullptr};
    std::vector<Fl_Box*> cell_widgets;
    int missing_value{0};

    Fl_Input *in{nullptr};
    Fl_Output *out{nullptr};
    Fl_Button *btn{nullptr};
    Fl_Button *btn_clear{nullptr};
    Fl_Grid *grid{nullptr};
    Fl_Box *ibox{nullptr};
    Fl_Box *title{nullptr};
    Fl_Box *info_box{nullptr};
    F8 f8;

  private:
    void create_cells();
    void clear_cells();
    static void button_callback (Fl_Widget *, void *data);
    static void button_clear(Fl_Widget *, void *data);
    void print_output();

  public:
    static constexpr const char *info = "ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г\n\
F8: В последовательности записаны целые числа от M до N ( M меньше N, M больше или равно 1) в произвольном порядке, \
но одно из чисел пропущено (остальные встречаются ровно по одному разу). N не превосходит 1000. \
Последовательность заканчивается числом 0. Определить пропущенное число.\n";

    WindowF8(int w, int h, const char *title = "ДЗ-9: F8");   
};