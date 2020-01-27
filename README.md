# ecc-rs

Experimental rust crate for performing basic elliptic curve operations, built on top of the
elliptic curve implementation in [*ring*](https://github.com/briansmith/ring).
The main aim of the project is to be able to provide an API for creating
elliptic curve points and manipulating them, including:

- affine/jacobian coordinates;
- point addition;
- point multiplication;
- point serialization/deserialization (both compressed and uncompressed);
- hash-to-curve (as specified in:
  [draft-irft-cfrg-hash-to-curve](https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-05)).
  
Currently, we only support the NIST P-256 and P-384 curves, as we are reliant on what is exposed
in *ring*.

Note that this is highly experimental code:

- it's **not** performant (as it converts between many data types a lot);
- it's **not** constant-time (as it uses BigUint and BigInt ops);
- it's actively **under-development** and there may be future API breaking changes;
- it has **not** undergone any sort of security review.

## Disclaimer

**This project is entirely experimental and has not undergone any sort of
security review. It SHOULD NOT be used in any production environment that seeks
to provide cryptographic security based on any of the primitives in this
library.**

## helpful commands

Build:
```
cargo build
```

Test:
```
cargo test
```

Docs:

```
cargo doc --lib --open
```
