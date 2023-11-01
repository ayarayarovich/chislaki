
def decimal_range(start, stop, increment):
    while start < stop: # and not math.isclose(start, stop): Py>3.5
        yield start
        start += increment

with open('for_two_dimensional_interpolation_sequential.txt', mode="w") as f:
  xx = []
  yy = []
  for x in decimal_range(-1, 1.0001, 1):
    xx.append(x)
    yy.append(x)
    print(x, file=f, end=' ')
  print(file=f)

  for y in yy:
    print(y, file=f, end=' ')
    for x in xx:
      print(f'{x**2 + y**2}', file=f, end=' ')
    print(file=f)




