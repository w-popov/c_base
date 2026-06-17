#include "gF1.hpp"
#include <FL/Fl_Button.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Multiline_Output.H>
#include <FL/Fl.H>
#include <string>
#include <vector>
#include <sstream>
#include "HW9.h"

MyWindow::MyWindow(int w, int h, const char *title)
    : Fl_Double_Window(w, h, title)
{
    begin();
    info = new Fl_Multiline_Output(50, 10, 450, 60);
    info->wrap(true);
    info->vertical_label_margin(5);
    info->value(INFO);
    info->textsize(12);
    info->color(FL_FREE_COLOR);
    info->box(FL_NO_BOX);

    in = new Fl_Input(50, 80, 450, 30);
    in->value("");
    in->textsize(14);
    
    out = new Fl_Output(50, 120, 450, 30);
    out->textsize(14);
    out->color(FL_LIGHT1);
    
    btn = new Fl_Button(50, 170, 150, 40, "Сортировать");
    btn->labelsize(14);
    btn->callback(btn_cb, this);
    end();
}

std::vector<int> &MyWindow::parse_input(const char *input)
{
    if (input == nullptr || *input == '\0')
    {
        array.clear();
        return array;
    }
    array.clear();
    std::stringstream ss(input);
    int number{0};
    while (ss >> number)
    {
        array.emplace_back(std::move(number));
    }
    return array;
}

void MyWindow::sort_arr(int *array, int size)
{
    sort_array(size, array);
}

void MyWindow::btn_cb(Fl_Widget*, void *data)
{
    MyWindow *win = static_cast<MyWindow *>(data);
    std::vector<int> parsed_array = win->parse_input(win->in->value());
    if (parsed_array.empty())
    {
        win->out->value("Введите числа для сортировки.");
        return;
    }
    win->sort_arr(parsed_array.data(), static_cast<int>(parsed_array.size()));
    std::string output;
    for (const auto &num : parsed_array)
    {
        output += std::to_string(num) + " ";
    }
    win->out->value(output.c_str());
}

int main ()
{
    Fl::set_font(FL_HELVETICA, "DejaVu Sans");
    MyWindow win(580, 300, "ДЗ-9: F1");
    win.size_range(580, 300);
    Fl::scheme("gtk+");
    win.show();
    return Fl::run();
}