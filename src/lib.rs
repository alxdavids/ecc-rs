//! The `ecc-rs` crate aims to provide an interface for performing low-level
//! elliptic curve operations. We use a modified version of the `ring` crate for
//! performing elliptic curve operations (in turn, using BoringSSL under the
//! hood). In terms of modifications to the *ring* library, it is has mostly been a
//! case of publicising some structs and functions that were previously private.
//!
//! The only thing that the library exposes currently are structs called
//! `AffinePoint` and `JacobianPoint` corresponding to affine and jacobian
//! respresentations of elliptic curve points, respectively. These structs are
//! used to implement necessary low-level functionality such as addition of
//! points and scalar multiplication. Point coordinates are currently stored as
//! BigUint values, but these values are converted to `ring` respresentations
//! before any operations are performed.
//!
//! While this library uses the well-travelled *ring* library, the layer that we
//! have added is not intended to be used in *anything* other than experimental
//! implementations. Moreover, this library is not developed with optimised
//! performance in mind and so anything that is performance critical may way want
//! to look elsewhere.

#![allow(
    missing_copy_implementations,
    missing_debug_implementations,
    non_camel_case_types,
    non_snake_case,
    unsafe_code
)]
#![deny(
    dead_code,
    missing_docs,
    unstable_features,
    unused_qualifications,
    variant_size_differences,
)]
#![forbid(
    anonymous_parameters,
    trivial_casts,
    trivial_numeric_casts,
    unused_extern_crates,
    unused_import_braces,
    unused_results,
)]
#![cfg_attr(feature = "internal_benches", allow(unstable_features), feature(test))]

pub mod point;