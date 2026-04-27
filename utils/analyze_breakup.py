import math

def analizar_energia_relativa(archivo_entrada):
    # Abrimos el archivo en modo lectura
    with open(archivo_entrada, 'r') as f:
        lineas = f.readlines()

    # 1. Extraer las masas (índices 2 y 3 corresponden a las líneas 3 y 4)
    m1 = float(lineas[2].split()[0])
    m2 = float(lineas[3].split()[0])

    # 2. Extraer parámetros de las mallas (líneas 7, 8 y 10)
    # Según la estructura de la línea: num_puntos, step_rad, min_deg, max_deg, step_deg
    malla_th1 = lineas[6].split()
    th1_min = float(malla_th1[2]) # Ángulo inicial en grados
    d_th1 = float(malla_th1[4])   # Paso en grados

    malla_th2 = lineas[7].split()
    th2_min = float(malla_th2[2])
    d_th2 = float(malla_th2[4])

    malla_phi = lineas[9].split()
    phi_min = float(malla_phi[2])
    d_phi = float(malla_phi[4])

    # Imprimimos una cabecera para presentar los resultados ordenadamente
    print(f"Masa Fragmento 1 (4He): {m1} u")
    print(f"Masa Fragmento 2 (3He): {m2} u")
    print(f"{'E1 (MeV)':<12} {'E2 (MeV)':<12} {'Sección Eficaz':<18} {'E_rel (MeV)':<12}")
    print("-" * 58)

    # 3. Procesar los datos (desde la línea 11 en adelante, que es el índice 10)
    for linea in lineas[10:]:
        datos = linea.split()
        
        # Saltamos líneas en blanco si existieran
        if not datos:
            continue
        
        # Extraemos los índices de las mallas (son enteros)
        i_th1 = int(datos[1])
        i_th2 = int(datos[2])
        i_phi = int(datos[3])
        
        # Manejamos la notación científica de Fortran ('D' en lugar de 'E')
        sigma_str = datos[4].replace('D', 'E')
        sigma = float(sigma_str)
        
        E1 = float(datos[5])
        E2 = float(datos[6])

        # 4. Condición para computar solo cuando la sección eficaz no es nula
        if sigma != 0.0:
            # Reconstruimos los ángulos en grados usando el índice (restamos 1 porque el índice empieza en 1)
            th1_deg = th1_min + (i_th1 - 1) * d_th1
            th2_deg = th2_min + (i_th2 - 1) * d_th2
            phi_deg = phi_min + (i_phi - 1) * d_phi

            # Convertimos a radianes para las funciones trigonométricas de Python
            th1 = math.radians(th1_deg)
            th2 = math.radians(th2_deg)
            phi = math.radians(phi_deg)

            # Calculamos el coseno del ángulo relativo entre ambos fragmentos (theta_12)
            cos_th12 = (math.sin(th1) * math.sin(th2) * math.cos(phi) + 
                        math.cos(th1) * math.cos(th2))

            # Calculamos la Energía Cinética Relativa
            term1 = m2 * E1
            term2 = m1 * E2
            term3 = 2 * math.sqrt(m1 * m2 * E1 * E2) * cos_th12
            
            E_rel = (term1 + term2 - term3) / (m1 + m2)

            # Mostramos el resultado en pantalla
            print(f"{E1:<12.6f} {E2:<12.6f} {sigma:<18.6e} {E_rel:<12.6f}")

# Ejecutamos la función asumiendo que el archivo se llama "fort.777"
analizar_energia_relativa('fort.777')
