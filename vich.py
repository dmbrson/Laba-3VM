from math import factorial
from numpy import cos, sin, pi, linspace

#вычисление производной от полинома Лагранжа первой степени в заданной точке.
def lagrange_derivative(num_point: int, points: list[tuple[float, float]], step: float) -> float:
    count_points = len(points) - 1
    result = 0
    for i, point in enumerate(points):
        point_value = point[1]

        #часть множителя разности
        def diff_mult_part(a: int, b: int) -> float:
            sub_mult = 1
            for j in range(a, b):
                sub_mult *= (i - j)
            return sub_mult

        diff_mult = diff_mult_part(0, i - 1 + 1) * diff_mult_part(i + 1, count_points + 1)

        #множитель сетки для каждой точки интерполяции
        def grid_mult_part(a: int, b: int) -> float:
            alt_mult = 0
            for j in range(a, b):
                sub_mult = 1
                for j1 in range(0, min(i, j) - 1 + 1):
                    sub_mult *= (num_point - j1)
                for j1 in range(min(i, j) + 1, max(i, j) - 1 + 1):
                    sub_mult *= (num_point - j1)
                for j1 in range(max(i, j) + 1, count_points + 1):
                    sub_mult *= (num_point - j1)
                alt_mult += sub_mult
            return alt_mult

        grid_mult = grid_mult_part(0, i - 1 + 1) + grid_mult_part(i + 1, count_points + 1)

        result += point_value / diff_mult * grid_mult

    return result / step


# Функция y = x - cos(x)
def target_function(x: float, deriv_order: int = 0):
    if deriv_order == 0:
        return x - cos(x)
    elif deriv_order == 1:
        return 1 + sin(x)
    return cos(x + ((deriv_order - 1) * pi) / 2)

#оценка теоретическоя погрешности
def theoretical_error_estimate(num_point: int, step: float,
                               count_points: int, range_limits: tuple[float, float], function) -> tuple[float, float]:
    pts = function(linspace(*range_limits, num=10 ** 3), count_points + 1)
    fct = factorial(count_points + 1)

    #определение минимального и максимального значения функции.
    result: tuple[float, float] = min(pts), max(pts)

    #вычисление погрешнонсти
    def calculate_sub_result(value: float) -> float:
        sub_res = 0
        for j in range(count_points + 1):
            sub_mult = 1
            for j1 in range(count_points + 1):
                if j1 != j:
                    sub_mult *= (num_point - j1)
            sub_res += sub_mult
        return value * sub_res * step ** count_points / fct
    return map(calculate_sub_result, result)


range_limits = (0.5, 1.0)
deriv_order = 1
count_points = 5
point_index = 5
step_size = (range_limits[1] - range_limits[0]) / count_points
#генерируем набор образцовых точек
sample_points = [(pt, target_function(pt)) for pt in linspace(*range_limits, count_points + 1)]
#значение производной в точке
lagrange_result = lagrange_derivative(point_index, sample_points, step_size)
#фактическо е значение
function_result = target_function(sample_points[point_index][0], deriv_order)
round = 10
#минимальная и максимальная погрешность
min_theoretical_error, max_theoretical_error = map(lambda num: round(num, round),
                                                   theoretical_error_estimate(point_index, step_size,
                                                                              count_points, range_limits,
                                                                              target_function))

output_file = "results.txt"
with open(output_file, "w") as file:
    file.write("Значение производной по методу Лагранжа:\t" + str(round(lagrange_result, round)) + "\n")
    file.write("Фактическое значение производной:\t" + str(round(function_result, round)) + "\n")
    file.write("Разница:\t" + str(round(abs(lagrange_result - function_result), round)) + "\n")
    file.write("\n")
    file.write("Минимальная теоретическая ошибка:\t" + str(min_theoretical_error) + "\n")
    file.write("Максимальная теоретическая ошибка:\t" + str(max_theoretical_error) + "\n")
    file.write("\n")
    file.write("Попадает ли ошибка в указанный интервал?:\t"
               + str(min_theoretical_error < abs(lagrange_result - function_result) <max_theoretical_error))