# Abstract Operators

```@meta
DocTestSetup = quote
	using ED_sectors
end
```

```@index
Pages=["abstract_operators.md"]
```


```@docs
TERM
*(::Number, ::TERM)
ABSTRACT_OP
*(::Number, ::ABSTRACT_OP)
+(::ABSTRACT_OP,::ABSTRACT_OP)
+(:: ABSTRACT_OP, :: TERM)
```