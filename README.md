# FastLink.jl
Fast Probabilistic Record Linkage for the Julia Language
## What is FastLink.jl

The purpose of FastLink.jl is to bring a fast record linkage package to the julia language. When attempting to match large datasets using existing libraries in R and Python, I found they can be very slow and succumb to issues with memory pressure. This implementation of the fastlink algorithm is intended to scale effeciently in parallel and be able to easily handle matches between tabular data that span millions of rows. 

[![Run tests](https://github.com/jw2249a/FastLink.jl/actions/workflows/test.yml/badge.svg)](https://github.com/jw2249a/FastLink.jl/actions/workflows/test.yml)
