---
anchor: intro
title: Introduction
order: 1
---

## Introduction

SSAdNDP [1] is an extension of the AdNDP method [2] to periodic systems, derived from the periodic implementation [3] of Natural Bond Orbital (NBO) analysis [4]. It allows interpretation of chemical bonding in systems with translational symmetry in terms of classical lone pairs and two-center bonds, as well as multi-center delocalized bonding. The bonding pattern is expressed as a set of *n*-center – 2-electron (*n*c-2e) bonds.

The original program was written by Alexander Boldyrev at Utah State University. The source code in this repository has been updated to replace deprecated routine-calling schemes and modernize the Makefile for compatibility with recent Intel compilers and MKL (tested with version 2022.1.1).

<div class="note" markdown="1">
**Citation required.** Anyone who downloads and uses the code must cite refs [1], [2], [3], and the NBO references [4].
</div>
