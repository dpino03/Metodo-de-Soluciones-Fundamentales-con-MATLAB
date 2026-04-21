# Método de las Soluciones Fundamentales (MFS) con MATLAB

Resolución numérica de problemas diferenciales mediante el **Método de las Soluciones Fundamentales** (Method of Fundamental Solutions, MFS), implementado en MATLAB. El repositorio contiene tanto el código como ejemplos reproducibles para cada tipo de problema abordado.

Este material es parte del Trabajo Fin de Grado *"Resolución de problemas diferenciales con el Método de Soluciones Fundamentales e implementación en MATLAB"*, defendido en la **Universidad de Sevilla** (Grado en Matemáticas, 2025) con calificación de **Sobresaliente (9,7/10)**.

---

## ¿Qué es el MFS?

El Método de las Soluciones Fundamentales es un método sin malla (*meshless*) para resolver problemas de contorno asociados a ecuaciones en derivadas parciales lineales. La idea central es aproximar la solución del problema como una combinación lineal de soluciones fundamentales del operador diferencial, centradas en puntos fuera del dominio (puntos fuente):

$$u(x) \approx \sum_{j=1}^{N} \alpha_j \, \Phi(x, y_j), \qquad y_j \notin \overline{\Omega}$$

Los coeficientes $\alpha_j$ se determinan imponiendo las condiciones de contorno en un conjunto de puntos de colocación sobre la frontera, lo que lleva a un sistema lineal.

Entre sus ventajas están la **ausencia de malla**, la **alta tasa de convergencia** para problemas con solución suave, y la **facilidad de implementación** frente a métodos como elementos finitos o diferencias finitas. Sus limitaciones principales son la sensibilidad a la colocación de los puntos fuente y el mal condicionamiento de la matriz del sistema, aspectos que se discuten en la memoria.

---

## Contenido del repositorio

El código está organizado por tipo de problema. Cada carpeta es autocontenida: incluye los scripts `.m` necesarios y ejemplos numéricos concretos.

| Carpeta | Problema | Operador |
|---|---|---|
| `Ecuación de Laplace` | Problemas de contorno estacionarios | $\Delta u = 0$ |
| `Ecuación de Ondas` | Propagación de ondas (problemas evolutivos) | $\partial_{tt} u - c^2 \Delta u = 0$ |
| `Ecuación del Calor` | Difusión (problemas parabólicos) | $\partial_t u - k \Delta u = 0$ |
| `Problema Elíptico` | Problemas elípticos más generales | $Lu = f$ |
| `Problema inverso Ecuación del Calor` | Identificación de condiciones desde datos | inverso asociado a la ecuación del calor |
| `Problemas inversos elípticos` | Recuperación de parámetros/fuentes | inverso asociado al problema elíptico |

Los problemas inversos son especialmente relevantes: son los casos donde el MFS muestra más ventajas frente a métodos clásicos, y también donde se manifiestan los retos de regularización que se analizan en el trabajo.

---

## Cómo ejecutar los ejemplos

**Requisitos:** MATLAB R2021a o superior (no se han usado toolboxes adicionales).

1. Clona el repositorio:
   ```bash
   git clone https://github.com/dpino03/Metodo-de-Soluciones-Fundamentales-con-MATLAB.git
   ```
2. Abre MATLAB y sitúate en la carpeta del problema que quieras ejecutar.
3. Ejecuta el script principal de esa carpeta. Los parámetros (número de puntos de colocación, distancia de los puntos fuente, etc.) están en las primeras líneas del script y se pueden modificar para explorar su efecto sobre la solución.

Las gráficas de la solución, el error y la convergencia se generan automáticamente al ejecutar los scripts.

---

## Estructura típica de cada carpeta

Dentro de cada carpeta se sigue un patrón similar:

- Script principal con el planteamiento del problema y los parámetros.
- Funciones auxiliares para construir la matriz del sistema, aplicar condiciones de contorno y calcular el error.
- Ejemplos concretos con solución analítica conocida, para poder validar numéricamente.

---

## Autor

**David Pino Molero**  
Graduado en Matemáticas por la Universidad de Sevilla.

- GitHub: [@dpino03](https://github.com/dpino03)
- LinkedIn: [davidpinomolero](https://linkedin.com/in/davidpinomolero)
- Email: dpinomolero@gmail.com
