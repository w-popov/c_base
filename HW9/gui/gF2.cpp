#include "gF2.hpp"
#include <FL/Fl_Button.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Grid.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Multiline_Output.H>
#include <FL/Fl.H>
#include <string>
#include <vector>
#include <sstream>
#include "HW9.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <FL/x.H>
#endif

WindowF2::WindowF2(int w, int h, const char *title)
    : Fl_Double_Window(w, h, title)
{
    grid = new Fl_Grid(0, 0, w, h);
    grid->layout(6, 6, 2, 2);

    info_box = new Fl_Box(0, 0, 550, 60, INFO2);
    info_box->align(FL_ALIGN_CENTER | FL_ALIGN_INSIDE | FL_ALIGN_WRAP);
    info_box->labelfont(FL_HELVETICA);
    info_box->labelsize(13);
    info_box->color(FL_FREE_COLOR);
    info_box->box(FL_NO_BOX);

    grid->widget(info_box, 0, 0, 1, 6, FL_GRID_HORIZONTAL | FL_GRID_FILL);
    grid->row_gap(0, 10);
    grid->row_height(0, 80);

    in = new Fl_Input(0, 0, 550, 30);
    in->textsize(14);
    grid->widget(in, 2, 2, 1, 3, FL_GRID_HORIZONTAL);
    grid->row_gap(2, 10);

    out = new Fl_Output(0, 0, 550, 30);
    out->color(FL_LIGHT1);
    out->textsize(14);
    grid->widget(out, 3, 2, 1, 3, FL_GRID_HORIZONTAL);
    grid->row_gap(3, 10);

    btn = new Fl_Button(0, 0, 200, 40, "Сортировать");
    btn->callback(button_callback, this);
    grid->widget(btn, 4, 2, 1, 1, FL_GRID_HORIZONTAL);
    grid->col_gap(2, 20);

    int cw[] = { 20,  0,  0, 20, 20, 20};
    int rw[] = { 20,  0,  0,  0,  0, 10};
    grid->col_weight(cw, 6);
    grid->row_weight(rw, 6);

    grid->end();    
    this->end();

    grid->color(fl_rgb_color(250, 250, 250));
    this->color(fl_rgb_color(250, 250, 250));

    this->resizable(grid);
    this->size_range(600, 250, 800, 300);
}

void WindowF2::button_callback(Fl_Widget*, void *data)
{
    WindowF2 *win = static_cast<WindowF2 *>(data);
    std::vector<int> &array = win->parse_input(win->in->value());
    if (array.empty())
    {
        win->out->value("Введите числа для сортировки.");
        return;
    }
    win->sort_eodd(array.data(), static_cast<int>(array.size()));
    win->print_output(array);
}

void WindowF2::print_output(const std::vector<int> &array)
{
    std::ostringstream oss;
    for (const auto &num : array)
    {
        oss << num << " ";
    }
    out->value(oss.str().c_str());
}

std::vector<int>& WindowF2::parse_input(const char *input)
{
    array.clear();
    std::istringstream iss(input);
    int number;
    while (iss >> number) {
        array.push_back(number);
    }
    return array;
}

/**
 * Вызов Си функции сортировки четных и нечетных чисел в массиве.
 */
void WindowF2::sort_eodd(int *array, int size)
{
    sort_even_odd(size, array);
}

int main(int argc, char **argv)
{
    Fl::scheme("gtk+");
    WindowF2 *win = new WindowF2(600, 300, "ДЗ-9: F2");
    win->show(argc, argv);
    
    Fl::check(); 

    #ifdef _WIN32
    HWND hwnd = fl_xid(*win);
    if (hwnd)
    {
        HICON hIcon = LoadIcon(GetModuleHandle(NULL), "MAINICON");
        if (hIcon)
        {
            SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)hIcon);
            SendMessage(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)hIcon);
        }
    }
    #endif
    return Fl::run();
}