#![allow(
    missing_copy_implementations,
    missing_debug_implementations,
    non_camel_case_types,
    non_snake_case,
    unsafe_code
)]
// `#[derive(...)]` uses `trivial_numeric_casts` and `unused_qualifications`
// internally.
#![deny(
    unstable_features, // Used by `internal_benches`
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
    warnings
)]
#![no_std]
#![cfg_attr(feature = "internal_benches", allow(unstable_features), feature(test))]

#[cfg(feature = "alloc")]
extern crate alloc;

pub mod bits;
pub mod c;
pub mod error;
pub mod limb;
pub mod rand;
pub mod io;
pub mod ec;

#[macro_use]
pub mod arithmetic;

#[macro_use]
pub mod debug;

mod sealed {
    /// Traits that are designed to only be implemented internally in *ring*.
    //
    // Usage:
    // ```
    // use crate::sealed;
    //
    // pub trait MyType: sealed::Sealed {
    //     // [...]
    // }
    //
    // impl sealed::Sealed for MyType {}
    // ```
    pub trait Sealed {}
}
