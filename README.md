# M-todo-SPPS-para-propagaci-n-de-ondas-en-medios-estratificados
Códigos empleados para el trabajo Aplicación del método de Series de Potencias del Parámetro Espectral para el 
análisis numérico de la propagación de ondas de radio en la ionosfera no homogénea estratificada sin aplicación de fuerzas a las cargas

awj.mat - Es una variable implementada para ser una matriz que contiene los números de modo de onda con sus frecuencias 83 < w < 95.987 rad/s respectivamente.

fastAlfasVg.m -  Es el código implementado para construir tanto la matriz awj.mat y obtener las velocidades de grupo de los modos.

ninteg.m - Código que permite realizar integrales numéricas que se usa para el proceso del método SPPS, puede usarse algún otro método 
para la resolución numérica de integrales, no forzosamente este, queda a gusto y/o decisión del usuario. Este debe ser agregado a los códigos que involucre al método SPPS.

metodoAnalitico.m - Código que resuelve la ecuación de dispersión por el método analítico con el modelo más simple de la atmósfera.

metodoSPPS.m - Código que resuelve la ecuación de dispersión con el método SPPS para el perfil estratificado de la ionosfera mostrado en el trabajo.

simpleMetodoSPPS.m - Código que resuelve la ecuación de dispersión por el método SPPS con el modelo más simple de la atmósfera.
