#pragma once

template <typename T>
int* fill_array (T& array)
{
    for (auto &item : array) {
        scanf("%d", &item);
    }
    return array.data();
}
