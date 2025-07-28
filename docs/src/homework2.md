# Dokumentacija za Domačo Nalogo 2

## Funkcionalnosti paketa

```@docs
de_casteljau
calculate_bezier_loop_area
```

---

## 1   Algoritem de Casteljau

[Algoritem de Casteljau](https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm) je stabilen in numerično robusten postopek za evalvacijo (in razbitje) Bézierjeve krivulje stopnje $n$.

Bézierova krivulja je podana parametriično s formulo:

$$\mathbf{B}(t)=\sum_{i=0}^{n}{n \choose i}(1-t)^{n-i}t^{i}\mathbf{P}_{i},\qquad t\in[0,1],$$

kjer so $\mathbf{P}_{i}$ kontrolne točke.

### 1.1   Rekurzivna interpolacija
Algoritem temelji na ponavljajoči se linearni interpolaciji med zaporednimi točkami:

$$\mathbf{P}^{(k)}_{i}(t)=(1-t)\,\mathbf{P}^{(k-1)}_{i}(t)\;+\;t\,\mathbf{P}^{(k-1)}_{i+1}(t),\quad k=1,\dots ,n,$$

pri čemer $\mathbf{P}^{(0)}_{i}=\mathbf{P}_{i}$.  
Po $n$ korakih ostane ena sama točka $\mathbf{P}^{(n)}_{0}(t)=\mathbf{B}(t)$.

### 1.2   Implementacija
V funkciji `de_casteljau(points, t)` se:

1. Kopira vhodni vektor koeficientov (1D projekcija),
2. V dveh vgrajenih zankah izvede zgornja interpolacija,
3. Vrne vrednost v $t$.

Zahvaljujoč zgolj seštevanju in množenju z realnimi števili je algoritem brez problemov z nestabilnostjo, ki pestijo direktno evalvacijo Bernstein-ove vsote za večje $n$.

---

## 2   Ploščina zaprtih Bézierjevih zank

### 2.1 Green's theorem
Naj bo $\mathbf{r}(t)=(x(t),y(t))$ gladka, zaprta krivulja z $t\in[0,1]$.  
[Green's theorem](https://en.wikipedia.org/wiki/Green%27s_theorem) poda ploščino

$$A=\frac12\oint_{\mathbf{r}}\bigl(x\,dy-y\,dx\bigr)=\frac12\int_{0}^{1}\!\bigl[x(t)\,y'(t)-x'(t)\,y(t)\bigr]\,dt.$$

Za Bézierjevo krivuljo dobimo polinom stopnje $n-1$. Integral torej lahko izračunamo eksaktno z [Gauss–Legendrejevo kvadraturo](https://dlmf.nist.gov/3.5). V nalogi uporabimo 20-točkovno kvadraturo, ki [omogoča zadostno natančnost](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature).

### 2.2   Koraki funkcije `calculate_bezier_loop_area`

1. Ločimo koordinate kontrolnih točk v vektorja `px`, `py`.
2. Izračunamo kontrolne točke odvodov $x'(t)$ in $y'(t)$ s formulo  
   $$\mathbf{P}'_{i}=n\bigl(\mathbf{P}_{i+1}-\mathbf{P}_{i}\bigr).$$
3. Za vsako vozlišče $\xi_{k}\in[-1,1]$ 20-točkovne kvadrature izvedemo preslikavo  
   $t=\tfrac12(\xi_{k}+1)$ in ovrednotimo
   $$f(t)=x(t)\,y'(t)-x'(t)\,y(t)$$
   z `de_casteljau`.
4. Seštejemo $w_{k}f(t_{k})$ (teža $w_{k}$) s faktorjem $\tfrac14$:  
   prvi $\tfrac12$ iz Green's theorema, drugi $\tfrac12$ zaradi preslikave intervala.

---

## 3   Primer uporabe

Kontrolni poligon (datoteka `task.jl`):

```julia
CONTROL_POLYGON = [
  ( 0.0, 0.0),
  ( 1.0, 1.0),
  ( 2.0, 3.0),
  ( 1.0, 4.0),
  ( 0.0, 4.0),
  (-1.0, 3.0),
  ( 0.0, 1.0),
  ( 1.0, 0.0)
]
area = calculate_bezier_loop_area(CONTROL_POLYGON) # approx: 1.966
```
