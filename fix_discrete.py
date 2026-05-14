import sys

with open("belam.f90", "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm)" in line:
        if i >= 400 and i <= 450:
            lines[i] = "           if (ikterm + ilv == lambda) then\n              term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) \n"
            lines[i+1] = "      &                 * ((-av)/(ac+av))**(ilv)\n           else\n              term_pref = 1.0d0\n           endif\n"

with open("belam.f90", "w") as f:
    f.writelines(lines)
print("Updated discrete term_pref")
