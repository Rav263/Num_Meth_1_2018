# Метод Гауса и прочее 
Первое задание 3-ого семестра 2018 года по Численным методам

# Методы
    1) Простой метод Гаусса
    2) Метод Гаусса с выбором главного элемента
    3) Поиск обратной матрицы системы
    4) Подсчёт определителя системы
    5) Подстчёт числа обусловленности матрицы системы
    6) Итерационный метод (метод верхней реаксации)
    7) Генерация матрицы по заданным законам

# Использование

Сборка ./make.sh

Дополнительные аргументы:

    1 -- Для подсчёта обратной матрицы
    2 -- Для вычисления Определителя
    3 -- Для простого метода Гаусса
    4 -- Для метода Гаусса с выбором главного элемента
    5 -- Для метода верхней релаксации
    6 -- Для вычисления числа обусловленности матрицы


Для ввода матрицы с клавиатуры используйте ./gauss 1 n args

    n -- размер матрицы
    args -- дополнительные аргументы

Для ввода матрицы из файла используйте ./gauss 2 file_name n args

    n -- размер матрицы
    args -- дополнительные аргументы
    file_name -- имя файла из которого вводится матрица

Для генерации матрицы по некоторым законам используйте ./gauss 3 x args

    x -- параметр функции генерации
    args -- дополнительные аргументы

# Законы генерации

    Aij = Q ^ (i + j) + 0.1 * (j - i), i != j
    Aij = (Q - 1) ^ (i + j), i == j

    fi = n * exp(x / i) * cos(x / i);

    Q = 0.995

