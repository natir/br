/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* local mod */
pub mod cli;
pub mod error;

pub mod correct;

pub fn build_methods<'a>(
    params: Option<Vec<String>>,
    solid: &'a pcon::solid::Solid,
    confirm: u8,
    max_search: u8,
) -> Vec<Box<dyn correct::Corrector + 'a>> {
    let mut methods: Vec<Box<dyn correct::Corrector + 'a>> = Vec::new();

    if let Some(ms) = params {
        for method in ms {
            match &method[..] {
                "one" => methods.push(Box::new(correct::One::new(solid, confirm))),
                "graph" => methods.push(Box::new(correct::Graph::new(&solid))),
                "greedy" => {
                    methods.push(Box::new(correct::Greedy::new(&solid, max_search, confirm)))
                }
                "gap_size" => methods.push(Box::new(correct::GapSize::new(&solid, confirm))),
                _ => unreachable!(),
            }
        }
    } else {
        methods.push(Box::new(correct::One::new(&solid, confirm)));
    }

    methods
}
