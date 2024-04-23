use crate::copy_to_bytes::*;

/* Struct for keeping track of bitmap header image information.

Reference:

Wikipedia - BMP file format
https://en.wikipedia.org/wiki/BMP_file_format
*/
pub struct BitmapInfo {
    pub total_file_size: i32, // In number of bytes
    pub data_offset: i32, // Offset to the image data, in number of bytes
    pub header_size: u32, // In number of bytes

    pub width: i32,  // Image dimensions in number of pixels
    pub height: i32, // Image dimensions, continued.

    pub plane_count: u16, // Just set this to 1
    pub bits_per_pixel: u16,
    pub compression_method: u32, // 0 = no compression
    pub image_size: u32, // Number of bytes used for image data
    pub horizontal_resolution: i32, // Horizontal image dpi
    pub vertical_resolution: i32, // Vertical image dpi
    pub color_palette_count: u32, // Use 16777216 colors
    pub important_colors_count: u32, // Set this to 0
}


pub fn fill_bitmap_header(
    bytes: &mut [u8], info: BitmapInfo) {
    bytes[0] = b'B';
    bytes[1] = b'M';
    let mut offset: usize = 2;
    copy_i32(bytes, info.total_file_size, &mut offset);
    offset = 10;
    copy_i32(bytes, info.data_offset, &mut offset);
    copy_u32(bytes, info.header_size, &mut offset);
    copy_i32(bytes, info.width, &mut offset);
    copy_i32(bytes, info.height, &mut offset);
    copy_u16(bytes, info.plane_count, &mut offset);
    copy_u16(bytes, info.bits_per_pixel, &mut offset);
    copy_u32(bytes, info.compression_method, &mut offset);
    copy_u32(bytes, info.image_size, &mut offset);
    copy_i32(bytes, info.horizontal_resolution, &mut offset);
    copy_i32(bytes, info.vertical_resolution, &mut offset);
    copy_u32(bytes, info.color_palette_count, &mut offset);
    copy_u32(bytes, info.important_colors_count, &mut offset);
    // println!("{}", offset);
}

/* Write the contents of the data array to a bitmap file. Based
on the examples given in the Rust documentation:
https://doc.rust-lang.org/std/fs/struct.File.html
*/
pub fn make_bitmap_file(filename: std::string::String,
                        data: &mut [u8]) -> std::io::Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(filename)?;
    file.write_all(&data)?;
    Ok(())
}
