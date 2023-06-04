/**********************************************************************************
 * @file wav_reader.h
 * @brief WAV file reader utility functions
 *
 * This header file provides utility functions for reading audio data from WAV files.
 * It includes a structure definition for the WAV file header and a function to read
 * the audio samples from a WAV file into a float array.
 *
 * The `read_wav` function reads the audio samples from a given WAV file and returns
 * the number of samples and a pointer to a dynamically allocated float array that
 * holds the real part of the samples. The imaginary part of the audio samples is
 * assumed to be 0, as the WAV file format primarily stores real audio data.
 *
 * `````````````````````````````````````````````````````````````````````````````````
 * @date 2023-06-04
***********************************************************************************/

use std::fs::File;
use std::io::Read;
use std::mem::size_of;
use std::slice;

#[repr(C)]
struct Header {
    riff: [u8; 4],
    overall_size: u32,
    wave: [u8; 4],
    fmt_chunk_marker: [u8; 4],
    length_of_fmt: u32,
    format_type: u16,
    channels: u16,
    sample_rate: u32,
    byterate: u32,
    block_align: u16,
    bits_per_sample: u16,
    data_chunk_header: [u8; 4],
    data_size: u32,
}

pub fn read_wav(file_name: &str) -> (u32, Vec<f32>) {
    let mut file = File::open(file_name).expect("Failed to open file");

    let mut header: Header = unsafe { std::mem::zeroed() };
    file.read_exact(unsafe {
        slice::from_raw_parts_mut(
            &mut header as *mut _ as *mut u8,
            size_of::<Header>(),
        )
    })
    .expect("Failed to read header");

    let num_samples = header.data_size / (header.channels * (header.bits_per_sample / 8)) as u32;
    let mut data: Vec<f32> = Vec::with_capacity(num_samples as usize);

    for _ in 0..num_samples {
        for _ in 0..header.channels {
            let mut data_in_channel: i16 = 0;
            file.read_exact(unsafe {
                slice::from_raw_parts_mut(
                    &mut data_in_channel as *mut _ as *mut u8,
                    size_of::<i16>(),
                )
            })
            .expect("Failed to read data_in_channel");
            data.push(data_in_channel as f32);
        }
    }

    (num_samples, data)
}
