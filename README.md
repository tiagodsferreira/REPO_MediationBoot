# Parte 1: Introdução à mediação e efeitos indiretos

Em estatística, um modelo de mediação permite identificar e explicar um **efeito da variável independente (X) sobre a variável dependente (Y) por meio da inclusão de uma terceira variável denominada variável mediadora (M)**. Segundo Baron e Kenny, uma variável é designada de mediadora "*na medida em que ela explica/é responsável pela relação entre o preditor e o resposta*" (Figura 1).  

Na figura 1 são denotados os efeitos intervenientes num modelo mediação em que *c* representa o efeito total, isto é, o efeito causado por X em Y, *cp* designa-se por efeito parcial, este considera a inclusão da variável mediadora no modelo, correspondendo à relação direta entre a variável independente e dependente, mantendo fixo M. Ou seja, é o efeito de X em Y devido a outras causas que não a mediadora. Pelo facto de considerar a mediadora, *cp* difere de *c*, pois correspondem a relações diferentes. O efeito indireto envolve as letras *a* e *b*. Enquanto *a* representa o efeito da variável independente na variável mediadora, *b* representa o efeito da variável mediadora na variável dependente, mantendo fixa a variável independente. O efeito total é a soma dos efeitos diretos e indiretos.

Neste modelo de mediação, há dois caminhos para a variável dependente. A variável independente deve prever a variável dependente, e a variável independente deve prever o mediador. A mediação é testada através de três regressões:

\begin{equation}
\tag{1}
Y = b1 + \textit{c}X + e1
\end{equation}


\begin{equation}
\tag{2}
M = b2 + \textit{a} X + e2 
\end{equation}

\begin{equation}
\tag{3}
Y= b3 + \textit{b}M + \textit{cp}X + e3 
\end{equation}

A metodologia mais utilizada para a análise de mediação é conhecida por métodos dos quatro passos. Genericamente analisam as equações de regressão estimadas, de forma a verificar, em cada passo, se os seus coeficientes são estatisticamente significativos. De uma forma sucinta:

* verificar a significância estatística do coeficiente *c*, através da equação (1);  
* verificar a significância estatística do coeficiente *a*, através da equação (2);  
* verificar a significância estatística do coeficiente *b*, através da equação (3), quando fixo X, .
* verificar a significância estatística do coeficiente *cp*, através da equação (3), quando fixo M.  

Dizemos que a mediação é total quando a estimativa de *cp* é nula, indicando que o efeito da variável independente na variável dependente é explicado apenad pela variável mediadora. Na mediação parcial a estimativa de *cp* é inferior à estimativa de *c* (em valor absoluto).  

Porém, outros autores como por exemplo Preacher e Hayes (2010), utilizaram métodos alternativos para estimar os coeficientes e verificar a sua significância. Para estimar o efeito de mediação (indireto) é referido o método do produto dos coeficientes (*ab*), acompanhado pelo respetivo teste de significância. Este teste de significância conjunta para o efeito indireto definido como o produto de coeficientes é baseado na distribuição do produto de coeficientes. No entanto, o método bootstrap é o mais recomendado.

A título de curiosidade, deixamos aqui alguns exemplos práticos que utilizam a análise de mediação:  
- https://www.sciencedirect.com/science/article/pii/S231472101500002X  
- https://www.thelancet.com/pdfs/journals/lanplh/PIIS2542-5196(21)00235-7.pdf  
- https://www.thelancet.com/pdfs/journals/lanpub/PIIS2468-2667(20)30292-9.pdf  
- https://www.thelancet.com/pdfs/journals/eclinm/PIIS2589-5370(21)00483-1.pdf  
- https://eng.uber.com/mediation-modeling/  
# 
