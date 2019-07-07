#!/bin/bash
ls ../src/model_restraint_*.py | while read file; do mod9.15 $file; done
