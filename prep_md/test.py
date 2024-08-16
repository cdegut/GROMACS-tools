import default_parameter

for param in default_parameter.md:
    print(f"{param}\t\t=\t{default_parameter.md[param][0]}\t;", end="")
    if len(default_parameter.md[param]) > 1:
        for i in range(1, len(default_parameter.md[param])):
            print(default_parameter.md[param][i], end="")
    print()