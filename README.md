This is a local ensemble Kalman filter implementation for [LightDA](https://github.com/LightDA-assim/lightda-core). It is based heavily on the implementation from [PDAF](http://pdaf.awi.de/trac/wiki).

Note that because the source code derives from PDAF, this component uses the [GNU LGPL v3](https://choosealicense.com/licenses/lgpl-3.0/) license rather than the [MIT](https://choosealicense.com/licenses/mit/) license used by LightDA itself. This places additional restrictions on its use, the primary one being that derivative works must release their source code under the same license.

## Compiling

The simplest way to build the lenkf-rsm component is using the superbuild, which downloads and compiles LightDA and its dependencies compiling the lenkf-rsm component. It is invoked using CMake as follows:

```bash
mkdir build
cd build
cmake ../superbuild
make
```

To build lightda-lenkf-rsm alone using previously installed copies of LightDA, replace the above ```cmake ../superbuild``` with ```cmake ..```.