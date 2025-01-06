# GDCelestial

Celestial physics implementation for Godot 4.

## Currently supported

- [x] 2D conic sections.
- [ ] 3D conic sections.
- [ ] 2D N-body.
- [ ] 3D N-body.

## Known problems

Current Kepler equation solver is a little dumb and do not converge in case of bad initial approximation.
Probably should stop reinventing the wheel and use a separate well-crafted library for orbital mechanics...

## Acknowledgments

Based on [gdextension template](https://github.com/nathanfranke/gdextension) by Aaron Franke and Nathan Franke.
Additional thanks to Bryan Weber for reference materials on Orbital Mechanics & Astrodynamics.
See [orbital-mechanics.space](https://orbital-mechanics.space/intro.html).
