(* -----------------------------------------------------------------

+------------------------------------------------------------------+
|                                                                  |
|                                                                  |
|                 ____ ___     _________ ___  _______              |
|                / __ `__ \   / ___/ __ `/ / / / ___/              |
|               / / / / / /  / /  / /_/ / /_/ (__  )               |
|              /_/ /_/ /_/  /_/   \__,_/\__, /____/                |
|                                      /____/                      |
|                              v 1.0                               |
|                                                                  |
+------------------------------------------------------------------+

This code was authored by David Lizarazo in 2022. It is a ray-optics
implementation for the design and analysis of single surface optics.

----------------------------------------------------------------- *)

BeginPackage["mrays`"]

VecNorm::usage =  "VecNorm[vector] return the cartesian norm of the given vector.";
VecNormalize::usage = "VecNormalize[vector] returns the given vector scaled to have unit length.";
LabeledArrow::usage = "LabeledArrow[origin, dest, label] returns an Arrow starting at origin and arriving at dest with a labeled position in its midpoint."

PhaseUnwrap = ResourceFunction["PhaseUnwrap"];

RandomInSphere::usage="RandomInSphere[sphereDensity, n, \"OutputAs\"-> {\"Cartesian\", \"Polar-Azimuth\"}] gives n points distributed in a sphere according to the density sphereDensity[\[Theta],\[Phi]]. 
sphereDensity is a function of angles \[Theta] and \[Phi], where \[Theta] is assumed to be the polar angle and \[Phi] the azimuth angle. sphereDensity does not need to be normalized. sphereDensity is interpreted in the sense that in a solid angle d\[CapitalOmega] = Sin[\[Theta]]d\[Theta] d\[Phi] the probability of having an event in that direction is sphereDensity[\[Theta],\[Phi]] d\[CapitalOmega].
OutputAs is an option that can be used to control the format of the output, if set to \"Polar-Azimuth\" then the returned list consists of 2-element lists with the first element being a polar angle and the second one an azimuth angle. If OutputAs is set to \"Cartesian\" then the returned list consists of 3-element lists which correspond to {x,y,z} coordinates on the unit sphere.";

Refractor2D::usage="Refractor2D[\[Sigma], \[Sigma]y, n1, n2, ys, {pz, py}, normalized] takes a curve \[Sigma][y] and its derivative \[Sigma]y[y], the abscissa ys that together with \[Sigma][ys] define a point on the surface \[Sigma][y], the refractive indices n1 and n2, and the origin of a ray {py,pz}. This function returns the refraction vector that describes the direction of refraction of a ray that originates in {py,pz} and encounters the surface at {ys,\[Sigma][ys]}.
For the provided value for ys the function returns a list with four elements {pVec, oVec, reVector, True||False} where pVec is the vector from which the ray originates, oVec is the vector at the surface where refraction occurs, refVector is a unit vector pointing in the direction of refraction of the ray on the surface, the last element is False if there was total internal reflection and True otherwise.
If the angle of incidence is such that there is total internal refraction, then the function returns the direction vector of that reflection instead.
This function makes not attempt at avoiding spurious intersections or more than one encounter of ray with the surface, and assumes that the medium with refractive index n1 is to the left of \[Sigma][y].
";

Refractor2DParametric::usage="Refractor2DParametric[\[Sigma]z, \[Sigma]y, \[Sigma]z\[Theta], \[Sigma]y\[Theta], \[Theta], n1, n2, {pz, py}, normalize] takes the parametrization of a curve {\[Sigma]z[\[Theta]], \[Sigma]y[\[Theta]]} and its derivatives {\[Sigma]z\[Theta][\[Theta]], \[Sigma]y\[Theta][\[Theta]]}, the vector {pz, py} from which a ray originates, and the refractive indices n1 and n2. 
For the provided value for \[Theta] (which defines a point in the curve) the function returns a list with four elements {pVec,oVec,reVector,True||False} where pVec is the vector from which the ray originates, oVec i the vector at the surface where refraction occurs, refVector is a unit vector pointing in the direction of refraction of the ray on the surface, the last element is False if there was total internal reflection and True otherwise.
If the angle of incidence is such that there is total internal refraction, then the function returns the direction vector of that reflection instead.
This function makes not attempt at avoiding spurious intersections or more than one encounter of ray with the surface, and assumes that the medium with refractive index n1 is to the left of \[Sigma][y].
";

Refractor3D::usage="Refractor3D[\[Sigma], \[Sigma]x, \[Sigma]y, n1, n2, {xs, ys}, {px, py, pz}, normalize] takes a surface \[Sigma][x, y], its partial derivatives \[Sigma]x[x, y] and \[Sigma]y[x, y], the coordinates {xs, ys}  of a point in the surface, the refractive indices n1 and n2, and the ray origin of a ray {px,py,pz}.
For the provided values {xs, ys} (which define a point in the surface together with \[Sigma][xs, ys]) the function returns a list with four elements {pVec,oVec,reVector,True||False} where pVec is the vector from which the ray originates, oVec i the vector at the surface where refraction occurs, refVector is a unit vector pointing in the direction of refraction of the ray on the surface, the last element is False if there was total internal reflection and True otherwise.
If the angle of incidence is such that there is total internal refraction, then the function returns that direction vector instead.";

Refractor3DParametric::usage="Refractor3DParametric[{{\[Sigma]x[\[Theta],\[Phi]],\[Sigma]y[\[Theta],\[Phi]],\[Sigma]z[\[Theta],\[Phi]]},{\[Sigma]x\[Theta][\[Theta],\[Phi]],\[Sigma]y\[Theta][\[Theta],\[Phi]],\[Sigma]z\[Theta][\[Theta],\[Phi]]},{\[Sigma]x\[Phi][\[Theta],\[Phi]],\[Sigma]y\[Phi][\[Theta],\[Phi]],\[Sigma]z\[Phi][\[Theta],\[Phi]]},{\[Theta],\[Phi]},n1,n2,{px,py,pz}] takes a surface parametrized by {\[Sigma]x\[Theta][\[Theta],\[Phi]],\[Sigma]y\[Theta][\[Theta],\[Phi]],\[Sigma]z\[Theta][\[Theta],\[Phi]]}, with derivatives {\[Sigma]x\[Theta][\[Theta],\[Phi]],\[Sigma]y\[Theta][\[Theta],\[Phi]],\[Sigma]z\[Theta][\[Theta],\[Phi]]} and {\[Sigma]x\[Phi][\[Theta],\[Phi]],\[Sigma]y\[Phi][\[Theta],\[Phi]],\[Sigma]z\[Phi][\[Theta],\[Phi]]}, refractive indices n1 and n2, values for {\[Theta],\[Phi]} (which specify a point in the surface), and a vector {px,py,pz} from which a light-ray is assumed to orginate from.
The function returns a list with four elements {pVec, oVec, reVector, True||False} where pVec is the vector from which the ray originates, oVec i the vector at the surface where refraction occurs, and refVector is a unit vector pointing in the direction of refraction of the ray on the surface.
If the angle of incidence is such that there is total internal refraction, then the function returns that direction vector instead.";

RayCollider2D::usage = "RayCollider2D[surfaceFun, yRange, {z0, y0}, numRays] calculates the points of intersection of a pencil of rays containing n rays that originates from {z0, y0} and intersects with the surface at point {zs, ys} as described by the function surfaceFun.
It returns a list whose elements are lists with two elements, the first element is the angle that the ray makes with the horizontal axis, the second elements is a list with three elements: the first element being zs, the second being ys, and the third element equal to the distance between {zs,ys} and {z0,y0}.
Function assumes that yRange is sorted.";

WithinInterval::usage = "WithinInterval[x, {xmin, xmax}] returns the value of the inequality (xmin <= x <= xmax).";

FindFirstRoot::usage = "FindFirstRoot[func, tmax, resolution] attempts to find the first root of func in the interval [0, tmax] by first detecting a change of sign and then finishing off root finding by using FindRoot. 
The size of the steps taken in the sign-changing part of the algorithm is equal to tmax/resolution.
If no change of sign is detected the function returns None.";

RayIntersection3D::usage = "RayIntersection3D[surfFun, {x0, y0, z0}, {\[Theta]ray, \[Phi]ray}, tmax, stepResolution] determines the intersection of a ray that originates at {x0, y0, z0} in direction {\[Theta]ray, \[Phi]ray} and the surface described by the function surfFun[x,y].
tmax is an estimate the dimension of a line with origin {x0, y0, z0} that should have definitely resulted in an itersection with the surface if there is any intersection at all, stepResolution is the number of steps that are attempted to find and intersection in the application of FindFirstRoot.";

RayCollider3D::usage = "RayCollider3D[surfaceFun, {x0, y0, z0}, n, {{xmin,xmax},{ymin,ymax}}, stepResolution, dzTol] calculates the intersection of n rays with origin {x0, y0, z0} that would intersect with the surface described by the function surfaceFun. {{xmin, xmax}, {ymin, ymax}} describes the extension the domain in the x-y plane. stepResolution is the number of steps taken in the implementation of FindFirstRoot and dzTol is the minimum discrepance in the intersections that is considered acceptable.";

Begin["`Private`"]

VecNorm[vec_] := (Sqrt[Total[vec^2]])
VecNormalize[vec_] := Simplify[(vec/VecNorm[vec])]
LabeledArrow[origin_, dest_, label_] := (
  arrow = Arrow[{origin, dest}];
  relVec = dest - origin;
  midpoint = (origin + dest)/2;
  \[Beta] = ArcTan[relVec[[1]], relVec[[2]]];
  \[Beta] = Mod[\[Beta], 2 \[Pi]];
  If[3 \[Pi]/2 > \[Beta] > \[Pi]/2, \[Beta] = \[Beta] - \[Pi]];
  text = 
   Text[label, midpoint, {0, -1}, {Cos[\[Beta]], Sin[\[Beta]]}];
  {arrow, text}
  )

Refractor2D[\[Sigma]_, \[Sigma]y_, n1_, n2_, ys_, {pz_, py_}, normalized_:False]:=(
  (*evaluate partial derivatives and other common factors in the final expression*)
  {\[Sigma]yv, \[Sigma]v} = {\[Sigma]y[ys], \[Sigma][ys]};
  (*This can be used to determine the condition for total internal reflection*)
  radicand = n2^2 - n1^2 (1.-(-pz+\[Sigma]v+py \[Sigma]yv-ys \[Sigma]yv)^2/(((py-ys)^2+(pz-\[Sigma]v)^2) (1+\[Sigma]yv^2)));
  pVec = {pz,py};
  oVec = {\[Sigma]v, ys};
  If[radicand<0.,
  ((* total internal reflection *)
    If[normalized,
      (factor = 1./(Sqrt[(py-ys)^2+(pz-\[Sigma]v)^2] (1.+\[Sigma]yv^2));
      reflectionVec = factor*{
        pz-2. py \[Sigma]yv+2. ys \[Sigma]yv-pz \[Sigma]yv^2+\[Sigma]v (-1.+\[Sigma]yv^2),
        -py+ys-2. (pz-\[Sigma]v) \[Sigma]yv+(py-ys) \[Sigma]yv^2
        };
      Return[{pVec, oVec, reflectionVec, False}];
      ),
      (reflectionVec = {
        pz-2. py \[Sigma]yv+2. ys \[Sigma]yv-pz \[Sigma]yv^2+\[Sigma]v (-1.+\[Sigma]yv^2),
        -py+ys-2. (pz-\[Sigma]v) \[Sigma]yv+(py-ys) \[Sigma]yv^2
        };
      Return[{pVec, oVec, reflectionVec, False}];
      )
    ]
  ),
  ((* refraction *)
    oneplusq = 1.+\[Sigma]yv^2;
    chunk1   = 1./n2*Sqrt[radicand/oneplusq];
    inormie  = 1./(Sqrt[(py-ys)^2+(pz-\[Sigma]v)^2] * oneplusq);
    nr=n1/n2;
    refractionVector={
      -nr( \[Sigma]yv (py-ys+(pz-\[Sigma]v) \[Sigma]yv))*inormie+chunk1,
      nr(-py+ys+(-pz+\[Sigma]v) \[Sigma]yv)*inormie-\[Sigma]yv chunk1
      };
    Return[{pVec, oVec, refractionVector, True}];
  )
  ]
);

Refractor2DParametric[\[Sigma]z_, \[Sigma]y_, \[Sigma]z\[Theta]_, \[Sigma]y\[Theta]_, \[Theta]_, n1_, n2_, {pz_,py_}, normalize_:False]:= ( \
  (*evaluate partial derivatives and other common factors in the final expression*) \
  {\[Sigma]zv, \[Sigma]yv, \[Sigma]z\[Theta]v, \[Sigma]y\[Theta]v} = {\[Sigma]z[\[Theta]], \[Sigma]y[\[Theta]], \[Sigma]z\[Theta][\[Theta]], \[Sigma]y\[Theta][\[Theta]]};
  (*This function value can be used to determine the condition for total internal reflection*) \
  radicand = n2^2-n1^2 (1-(pz \[Sigma]y\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]zv+(-py+\[Sigma]yv) \[Sigma]z\[Theta]v)^2/(((py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2) (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)));
  pVec = {pz, py};
  oVec = {\[Sigma]zv, \[Sigma]yv};

  If[radicand<0.,
    ((* total internal reflection *) \
    If[normalize, (
      factor = 1./(Sqrt[(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2));
      reflectionVec = factor*{-\[Sigma]y\[Theta]v^2 \[Sigma]zv+2 (-py+\[Sigma]yv) \[Sigma]y\[Theta]v \[Sigma]z\[Theta]v+\[Sigma]zv \[Sigma]z\[Theta]v^2+pz (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2),2 \[Sigma]y\[Theta]v (-pz+\[Sigma]zv) \[Sigma]z\[Theta]v-py (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2)+\[Sigma]yv (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2)};
      Return[{pVec, oVec, reflectionVec, False}];
      ),
      (
      reflectionVec = {
        -\[Sigma]y\[Theta]v^2 \[Sigma]zv+2 (-py+\[Sigma]yv) \[Sigma]y\[Theta]v \[Sigma]z\[Theta]v+\[Sigma]zv \[Sigma]z\[Theta]v^2+pz (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2),
        2 \[Sigma]y\[Theta]v (-pz+\[Sigma]zv) \[Sigma]z\[Theta]v-py (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2)+\[Sigma]yv (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2)
        };
      Return[{pVec, oVec, reflectionVec, False}];
      )
    ];
    ),
    ((* refraction *)
      root   = 1./n2*Sqrt[radicand/(\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)];
      factor = n1/n2 1./(Sqrt[(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2));
      refractionVector = {-((n1 \[Sigma]z\[Theta]v (py \[Sigma]y\[Theta]v-\[Sigma]yv \[Sigma]y\[Theta]v+(pz-\[Sigma]zv) \[Sigma]z\[Theta]v))/(n2 \[Kappa]^2 Sqrt[(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)))+(\[Kappa] \[Sigma]y\[Theta]v  root)/Sqrt[\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2],
      -((n1 \[Sigma]y\[Theta]v (py \[Sigma]y\[Theta]v-\[Sigma]yv \[Sigma]y\[Theta]v+(pz-\[Sigma]zv) \[Sigma]z\[Theta]v))/(n2 \[Kappa]^2 Sqrt[(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)))-(\[Kappa] \[Sigma]z\[Theta]v  root)/Sqrt[\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2]};
      Return[{pVec, oVec, refractionVector, True}];
    )
  ]
);

Refractor3D[\[Sigma]_, \[Sigma]x_, \[Sigma]y_, n1_, n2_, {xs_, ys_}, {px_, py_, pz_}, normalize_:False]:=(
  (*evaluate partial derivatives and other common factors in the final expression*)
  {\[Sigma]xv, \[Sigma]yv, \[Sigma]v} = {\[Sigma]x[xs,ys], \[Sigma]y[xs,ys], \[Sigma][xs,ys]};
  (*This can be used to determine the condition for total internal reflection*)
  radicand = (n2^2-n1^2 (1-(-pz+\[Sigma]v+(px-xs) \[Sigma]xv+py \[Sigma]yv-ys \[Sigma]yv)^2/(((px-xs)^2+(py-ys)^2+(pz-\[Sigma]v)^2) (1+\[Sigma]xv^2+\[Sigma]yv^2))));
  pVec = {px, py, pz};
  oVec = {xs, ys, \[Sigma]v};

  If[radicand<0,
    (
    eNorm = Sqrt[(px-xs)^2+(py-ys)^2+(pz-\[Sigma]v)^2];
    combo = 1./(eNorm (1+\[Sigma]xv^2+\[Sigma]yv^2));
    reflectionVec={
      combo((px-xs) \[Sigma]xv^2-2 \[Sigma]xv (pz-\[Sigma]v+(-py+ys) \[Sigma]yv)-(px-xs) (1+\[Sigma]yv^2)),
      combo(-py+ys+(-py+ys) \[Sigma]xv^2-2 (pz-\[Sigma]v) \[Sigma]yv+2 (px-xs) \[Sigma]xv \[Sigma]yv+(py-ys) \[Sigma]yv^2),
      combo(pz-2 (px-xs) \[Sigma]xv-pz \[Sigma]xv^2-2 py \[Sigma]yv+2 ys \[Sigma]yv-pz \[Sigma]yv^2+\[Sigma]v (-1+\[Sigma]xv^2+\[Sigma]yv^2))
      };
    Return[{pVec, oVec, reflectionVec, False}];
    ),
    (
    eNorm   = Sqrt[(px-xs)^2+(py-ys)^2+(pz-\[Sigma]v)^2];
    radical = 1./n2 Sqrt[radicand];
    combo1  = radical/Sqrt[1+\[Sigma]xv^2+\[Sigma]yv^2];
    combo2  = n1/n2  1/(eNorm (1+\[Sigma]xv^2+\[Sigma]yv^2));
    refractionVector = {
      -combo2 (\[Sigma]xv (pz-\[Sigma]v+(-py+ys) \[Sigma]yv)+(px-xs) (1+\[Sigma]yv^2))-\[Sigma]xv combo1,
       combo2 (-py+ys+(-py+ys) \[Sigma]xv^2+(-pz+\[Sigma]v) \[Sigma]yv+(px-xs) \[Sigma]xv \[Sigma]yv)-\[Sigma]yv combo1,
      -combo2 ((px-xs) \[Sigma]xv+(pz-\[Sigma]v) \[Sigma]xv^2+\[Sigma]yv (py-ys+(pz-\[Sigma]v) \[Sigma]yv))+combo1};
    Return[{pVec, oVec, refractionVector, True}];
    )
  ]
)

Refractor3DParametric[{\[Sigma]x_,\[Sigma]y_,\[Sigma]z_},{\[Sigma]x\[Theta]_,\[Sigma]y\[Theta]_,\[Sigma]z\[Theta]_},{\[Sigma]x\[Phi]_,\[Sigma]y\[Phi]_,\[Sigma]z\[Phi]_},{\[Theta]_,\[Phi]_},n1_,n2_,{px_,py_,pz_}]:=(
  (*evaluate partial derivatives and other common factors in the final expression*)
  {\[Sigma]xv, \[Sigma]yv, \[Sigma]zv}                         = {\[Sigma]x[\[Theta],\[Phi]],\[Sigma]y[\[Theta],\[Phi]],\[Sigma]z[\[Theta],\[Phi]]};
  {\[Sigma]x\[Theta]v, \[Sigma]y\[Theta]v, \[Sigma]z\[Theta]v} = {\[Sigma]x\[Theta][\[Theta],\[Phi]],\[Sigma]y\[Theta][\[Theta],\[Phi]],\[Sigma]z\[Theta][\[Theta],\[Phi]]};
  {\[Sigma]x\[Phi]v, \[Sigma]y\[Phi]v, \[Sigma]z\[Phi]v}       = {\[Sigma]x\[Phi][\[Theta],\[Phi]],\[Sigma]y\[Phi][\[Theta],\[Phi]],\[Sigma]z\[Phi][\[Theta],\[Phi]]};
  (* If this is negative, then total internal reflection occurs. *)
  radicand = n2^2-n1^2 (1-(pz \[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-pz \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v-\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v \[Sigma]zv+\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]zv-py \[Sigma]x\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]x\[Phi]v \[Sigma]yv \[Sigma]z\[Theta]v+px \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]xv \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+py \[Sigma]x\[Theta]v \[Sigma]z\[Phi]v-\[Sigma]x\[Theta]v \[Sigma]yv \[Sigma]z\[Phi]v-px \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v+\[Sigma]xv \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2/(((px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2) ((\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v)^2+(\[Sigma]x\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]z\[Phi]v)^2+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2)));
  pVec = {px, py, pz};
  oVec = {\[Sigma]xv, \[Sigma]yv, \[Sigma]zv};
  If[radicand<0 ,
    ( (* Total internal reflection *)
    reflectionVec = {(-2 (pz \[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-pz \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v-\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v \[Sigma]zv+\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]zv-py \[Sigma]x\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]x\[Phi]v \[Sigma]yv \[Sigma]z\[Theta]v+py \[Sigma]x\[Theta]v \[Sigma]z\[Phi]v-\[Sigma]x\[Theta]v \[Sigma]yv \[Sigma]z\[Phi]v) (-\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)-px (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)-(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))+\[Sigma]xv (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)-(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2)))/(Sqrt[(px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))),
    (\[Sigma]x\[Theta]v^2 \[Sigma]yv \[Sigma]y\[Phi]v^2+\[Sigma]yv \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v^2+\[Sigma]x\[Phi]v^2 (2 \[Sigma]y\[Theta]v (-pz+\[Sigma]zv) \[Sigma]z\[Theta]v+\[Sigma]yv (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2))-2 pz \[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v+2 \[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v \[Sigma]zv \[Sigma]z\[Phi]v+2 px \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v-2 \[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v-2 \[Sigma]yv \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v-\[Sigma]x\[Theta]v^2 \[Sigma]yv \[Sigma]z\[Phi]v^2-2 px \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v^2+2 \[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v^2+\[Sigma]yv \[Sigma]y\[Theta]v^2 \[Sigma]z\[Phi]v^2-2 \[Sigma]x\[Phi]v ((px-\[Sigma]xv) \[Sigma]z\[Theta]v (\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)-\[Sigma]x\[Theta]v (pz-\[Sigma]zv) (\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v \[Sigma]yv (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v-\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v))-py (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2+\[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (-2 \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+2 \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2-\[Sigma]z\[Phi]v^2)))/(Sqrt[(px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))),
    (-\[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v^2 \[Sigma]zv-2 px \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v+2 \[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v+\[Sigma]y\[Phi]v^2 \[Sigma]zv \[Sigma]z\[Theta]v^2+\[Sigma]x\[Phi]v^2 (-\[Sigma]y\[Theta]v^2 \[Sigma]zv+2 (-py+\[Sigma]yv) \[Sigma]y\[Theta]v \[Sigma]z\[Theta]v+\[Sigma]zv \[Sigma]z\[Theta]v^2)-2 py \[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v+2 \[Sigma]x\[Theta]v^2 \[Sigma]yv \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v+2 px \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v-2 \[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v-2 \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]zv \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v+\[Sigma]x\[Theta]v^2 \[Sigma]zv \[Sigma]z\[Phi]v^2+\[Sigma]y\[Theta]v^2 \[Sigma]zv \[Sigma]z\[Phi]v^2+pz (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2-\[Sigma]z\[Theta]v^2)-(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2+\[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (-2 \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+2 \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2-\[Sigma]z\[Phi]v^2))+2 \[Sigma]x\[Phi]v ((px-\[Sigma]xv) \[Sigma]y\[Theta]v (\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]zv+(py-\[Sigma]yv) \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v+\[Sigma]z\[Theta]v (py \[Sigma]y\[Phi]v-\[Sigma]yv \[Sigma]y\[Phi]v-\[Sigma]zv \[Sigma]z\[Phi]v))))/(Sqrt[(px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2)))};
    Return[{pVec,oVec,reflectionVec,False}];
    ),
    ( (* Ordinary refraction *)
    radical = Sqrt[radicand]/n2;
    refractionVector = {-((n1 ((pz \[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-pz \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v-\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v \[Sigma]zv+\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]zv-py \[Sigma]x\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]x\[Phi]v \[Sigma]yv \[Sigma]z\[Theta]v+py \[Sigma]x\[Theta]v \[Sigma]z\[Phi]v-\[Sigma]x\[Theta]v \[Sigma]yv \[Sigma]z\[Phi]v) (-\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)+px (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))-\[Sigma]xv (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))))/(n2 Sqrt[(px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))))+((-\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v) radical)/Sqrt[(\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v)^2+(\[Sigma]x\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]z\[Phi]v)^2+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2],(n1 (\[Sigma]x\[Theta]v^2 \[Sigma]yv \[Sigma]y\[Phi]v^2+\[Sigma]yv \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v^2+\[Sigma]x\[Phi]v^2 \[Sigma]y\[Theta]v (\[Sigma]yv \[Sigma]y\[Theta]v+(-pz+\[Sigma]zv) \[Sigma]z\[Theta]v)-pz \[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v+\[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v \[Sigma]zv \[Sigma]z\[Phi]v+px \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v-\[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v-2 \[Sigma]yv \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v-px \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v^2+\[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]z\[Phi]v^2+\[Sigma]yv \[Sigma]y\[Theta]v^2 \[Sigma]z\[Phi]v^2-py (\[Sigma]x\[Phi]v^2 \[Sigma]y\[Theta]v^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]x\[Theta]v^2 \[Sigma]y\[Phi]v^2+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2)+\[Sigma]x\[Phi]v (-2 \[Sigma]x\[Theta]v \[Sigma]yv \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v-(px-\[Sigma]xv) \[Sigma]z\[Theta]v (\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v (pz-\[Sigma]zv) (\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v))))/(n2 Sqrt[(px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2)))+((\[Sigma]x\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]z\[Phi]v) radical)/Sqrt[(\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v)^2+(\[Sigma]x\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]z\[Phi]v)^2+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2],-((n1 (px \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v-\[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v+pz \[Sigma]y\[Phi]v^2 \[Sigma]z\[Theta]v^2-\[Sigma]y\[Phi]v^2 \[Sigma]zv \[Sigma]z\[Theta]v^2+\[Sigma]x\[Phi]v^2 \[Sigma]z\[Theta]v (-\[Sigma]yv \[Sigma]y\[Theta]v+(pz-\[Sigma]zv) \[Sigma]z\[Theta]v)-\[Sigma]x\[Theta]v^2 \[Sigma]yv \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v-px \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v+\[Sigma]xv \[Sigma]x\[Theta]v \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Phi]v-2 pz \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v+2 \[Sigma]y\[Theta]v \[Sigma]y\[Phi]v \[Sigma]zv \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v+pz \[Sigma]x\[Theta]v^2 \[Sigma]z\[Phi]v^2+pz \[Sigma]y\[Theta]v^2 \[Sigma]z\[Phi]v^2-\[Sigma]x\[Theta]v^2 \[Sigma]zv \[Sigma]z\[Phi]v^2-\[Sigma]y\[Theta]v^2 \[Sigma]zv \[Sigma]z\[Phi]v^2+py (\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v) (\[Sigma]x\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Phi]v (2 \[Sigma]x\[Theta]v (-pz+\[Sigma]zv) \[Sigma]z\[Theta]v \[Sigma]z\[Phi]v+(px-\[Sigma]xv) \[Sigma]y\[Theta]v (-\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v \[Sigma]yv (\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v+\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v))))/(n2 Sqrt[(px-\[Sigma]xv)^2+(py-\[Sigma]yv)^2+(pz-\[Sigma]zv)^2] (\[Sigma]x\[Phi]v^2 (\[Sigma]y\[Theta]v^2+\[Sigma]z\[Theta]v^2)+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2-2 \[Sigma]x\[Theta]v \[Sigma]x\[Phi]v (\[Sigma]y\[Theta]v \[Sigma]y\[Phi]v+\[Sigma]z\[Theta]v \[Sigma]z\[Phi]v)+\[Sigma]x\[Theta]v^2 (\[Sigma]y\[Phi]v^2+\[Sigma]z\[Phi]v^2))))+((-\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v+\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v) radical)/Sqrt[(\[Sigma]x\[Phi]v \[Sigma]y\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]y\[Phi]v)^2+(\[Sigma]x\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]x\[Theta]v \[Sigma]z\[Phi]v)^2+(\[Sigma]y\[Phi]v \[Sigma]z\[Theta]v-\[Sigma]y\[Theta]v \[Sigma]z\[Phi]v)^2]};
    Return[{pVec,oVec,refractionVector,True}];
    )
  ]
  )

Options[RandomInSphere]={"OutputAs"->"Polar-Azimuth"};
RandomInSphere[sphereDensity_, numPoints_, OptionsPattern[]] := Module[
  {maxP, toRando, \[Theta], \[Phi], rando\[Theta], rando\[Phi], randoTops, randoFlips, goodOnes, randos, randoPoints, outputType},
  (
    maxP = NMaximize[{sphereDensity[\[Theta],\[Phi]]Sin[\[Theta]],0<=\[Theta]<=\[Pi],0<\[Phi]<=2\[Pi]},{\[Theta],\[Phi]}][[1]];
    randoPoints = {};
    outputType = OptionValue["OutputAs"];
    While[Length[randoPoints] < numPoints,
      (
      toRando       = numPoints - Length[randoPoints];
      rando\[Theta] = RandomReal[{0, \[Pi]},toRando];
      rando\[Phi]   = RandomReal[{0, 2\[Pi]},toRando];
      randoTops     = MapThread[sphereDensity[#1,#2] * Sin[#1]&, {rando\[Theta],rando\[Phi]}];
      randoFlips    = RandomReal[{0,maxP},toRando];
      goodOnes      = MapThread[Less,{randoFlips,randoTops}];
      randos        = Pick[Transpose[{rando\[Theta],rando\[Phi]}],goodOnes];
      randoPoints   = Join[randoPoints,randos];
      )
    ];
    If[Length[randoPoints]>numPoints,
      randoPoints=randoPoints[[;;numPoints]];
    ];
    (*Return according to output type.*)
    Switch[
      outputType,
      "Polar-Azimuth",
      Return[randoPoints],
      "Cartesian",
      Return[{Sin[#[[1]]]Cos[#[[2]]],Sin[#[[1]]]Sin[#[[2]]],Cos[#[[1]]]}&/@randoPoints]
    ];
)
]

WithinInterval[x_, {xmin_, xmax_}] := (Return[xmin <= x <= xmax]);

(* Options[RayCollider2D,{"Method" -> "Piecewise Interpolation"}] *)

RayCollider2D[surfaceFun_, yRange_, {z0_,y0_}, numRays_, OptionsPattern[]]:=(
  Off[Interpolation::inhr];
  {ymin, ymax} = {yRange[[1]],yRange[[-1]]};
  zs = {surfaceFun[#],#}&/@yRange;
  (* Calculate the angles formed wrt the point {z0,y0} according to an even sampling on the y-axis *)
  \[Beta]s = ArcTan[#[[1]]-z0,#[[2]]-y0]&/@zs;
  \[Beta]s = PhaseUnwrap[\[Beta]s];
  (* Split the ordinates in chunks that are monotonic. *)
  mono\[Beta]s = SequenceCases[\[Beta]s,x_/;Less@@x||Greater@@x];
  (* Collect the correspondig values of y. *)
  splitLengths = Length/@mono\[Beta]s;
  monoys = MapThread[Take[yRange,{#1,#2}]&,{Prepend[Drop[Accumulate[splitLengths]+1,-1],1],Accumulate[splitLengths]}];

  (* Create an evenly spaced range for the angles of rays in pencil *)
  {\[Beta]min,\[Beta]max}=MinMax[\[Beta]s];
  \[Beta]Range = Range[\[Beta]min,\[Beta]max,(\[Beta]max-\[Beta]min)/(numRays-1)];
  (* Create a series of interpolations for each of the chunks *)
  (* These interpolations take in launching angles from center of pencil and give values of y that correspond to points in the intersecting curve *)
  interpols = MapThread[Interpolation[Transpose[{#1,#2}]]&,{mono\[Beta]s,monoys}];
  (* Evaluate the interpolations over the values of \[Beta], if the value of \[Beta] is outside the domain of the interpolation give None to filter after*)
  collisions = Table[{\[Beta],If[WithinInterval[\[Beta],#[[1]][[1]]],#[\[Beta]], None]&/@interpols},{\[Beta],\[Beta]Range}];
  (* Filter on None output *)
  collisions = Table[{row[[1]],Select[row[[2]],Not[#===None]&]},{row,collisions}];
  (* Evaluate the ordinates corresponding to the values of y *)
  collisions = {#[[1]],Transpose[{surfaceFun/@#[[2]],#[[2]]}]}&/@collisions;
  (* Calculate the distance to {z0,y0}, this allows selecting the first intersection *)
  collisions = Table[{row[[1]],{#[[1]],#[[2]],Norm[#-{z0,y0}]}&/@row[[2]]},{row,collisions}];
  (* Leave out empty results ? *)
  collisions = Select[collisions,Length[#[[2]]]!=0&];
  (* Sort according to distance and choose the first one *)
  collisions = Table[{row[[1]],SortBy[row[[2]],Last][[1]]},{row,collisions}];
  Return[collisions]
)

FindFirstRoot[fun_, tmax_, resolution_] := (
  Module[{t},
   (
    dt = tmax/resolution;
    tee = 0;
    startSign = Sign[fun[tee]];
    sign = startSign;
    While[sign == startSign && tee < tmax,
     (tee = tee + dt;
      sign = Sign[fun[tee]];
      )];
    froot = If[sign != startSign,
      var /. FindRoot[Evaluate[fun[var]] == 0, {var, tee}],
      Null
      ];
    Return[froot];
    )
   ]
  )

RayIntersection3D[surfFun_, {x0_, y0_, z0_}, {\[Theta]ray_, \[Phi]ray_}, tmax_, stepResolution_] := (
  Module[{t, eqn},
   (
    {xray, yray, zray} = {
        Sin[\[Theta]ray]*Cos[\[Phi]ray], 
        Sin[\[Theta]ray]*Sin[\[Phi]ray],
        Cos[\[Theta]ray]} // N;
    eqn[var_] := Evaluate[(surfFun[x0 + var * xray, y0 + var*yray] - (z0 + var*zray))];
    tsol = FindFirstRoot[eqn, tmax, stepResolution];
    bestVDist = surfFun[x0 + tsol* xray, y0 + tsol*yray] - (z0 + tsol*zray);
    Return[{tsol, tsol*{xray, yray, zray} + {x0, y0, z0}, Abs[bestVDist]}]
    )
   ]
  )

RayCollider3D[surfFun_, {x0_, y0_, z0_}, numRays_, {{xmin_, xmax_}, {ymin_, ymax_}}, stepResolution_, dzTol_ : 0.01] := (
  Module[{rayPencil, collisions, remainingRays, conversionFraction, extraRays, moCollisions},
   (
    (* Generate a set of ray-angles using a uniform distribution *)
    rayPencil = RandomInSphere[1 &, numRays];
    (*Take a quick look at the surface*)
    surfaceSample = Outer[surfFun, Range[xmin, xmax, (xmax - xmin)/100], Range[ymin, ymax, (ymax - ymin)/100]];
    surfaceSample = Flatten[surfaceSample];
    (* Use that to determine the extent of solution domain *)
    {zmin, zmax} = MinMax[surfaceSample];
    (* Extend to include the ray origin. *)
    zmin = Min[zmin, z0]; zmax = Max[zmax, z0];
    (* How big a ray could fit. *)
    tmax = Sqrt[(xmax - xmin)^2 + (ymax - ymin)^2 + (zmax - zmin)^2];

    collisions = RayIntersection3D[surfFun, {x0, y0, z0}, #, tmax, stepResolution] & /@ rayPencil;
    collisions = Select[collisions, Not[#[[1]] === Null] &];
    collisions = Select[collisions, ((#[[1]] > 0) && (xmin <= #[[2]][[1]] <= xmax) && (ymin <= #[[2]][[2]] <= ymax)) && (#[[3]] <= dzTol ) &];
    (* Estimate how many rays actually reached the surface *)
    remainingRays = numRays - Length[collisions];
    (* Estimate the fraction that was successful*)
    conversionFraction = Length[collisions]/numRays;
    (* Iterate until the required number of rays has been collected *)
    While[Length[collisions] < numRays,
     (remainingRays = numRays - Length[collisions];
      extraRays = Round[remainingRays/conversionFraction];
      rayPencil = RandomInSphere[1 &, extraRays];
      moCollisions = RayIntersection3D[surfFun, {x0, y0, z0}, #, tmax, stepResolution] & /@ rayPencil;
      moCollisions = Select[moCollisions, ((#[[1]] > 0 || True) && (xmin <= #[[2]][[1]] <= xmax) && (ymin <= #[[2]][[2]] <= ymax)) && (#[[3]] <= dzTol ) &];
      collisions = Join[collisions, moCollisions];)
     ];
    collisions = collisions[[;; numRays]];
    Return[collisions];
    )
   ]
  )

End[]

EndPackage[]