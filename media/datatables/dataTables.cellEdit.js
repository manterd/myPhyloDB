/*! CellEdit 1.0.0
 * ©2016 Elliott Beaty - datatables.net/license
 */

/**
 * @summary     CellEdit
 * @description Make a cell editable when clicked upon
 * @version     1.0.0
 * @file        dataTables.editCell.js
 * @author      Elliott Beaty
 * @contact     elliott@elliottbeaty.com
 * @copyright   Copyright 2016 Elliott Beaty
 *
 * This source file is free software, available under the following license:
 *   MIT license - http://datatables.net/license/mit
 *
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 *
 * For details please refer to: http://www.datatables.net
 */

jQuery.fn.dataTable.Api.register('MakeCellsEditable()', function (settings) {
    var table = this.table();

    jQuery.fn.extend({

        // UPDATE
        updateEditableCell: function (callingElement) {
            var row = table.row($(callingElement).parents('tr'));
            var cell = table.cell($(callingElement).parent());
            var columnIndex = cell.index().column;
            var inputField =getInputField(callingElement);

            // Update
            var newValue = inputField.val();
            if (!newValue && ((settings.allowNulls) && settings.allowNulls != true)) {
                // If columns specified
                if (settings.allowNulls.columns) {
                    // If current column allows nulls
                    if (settings.allowNulls.columns.indexOf(columnIndex) > -1) {
                        _update(newValue);
                    } else {
                        _addValidationCss();
                    }
                    // No columns allow null
                } else if (!newValue) {
                    _addValidationCss();
                }
                //All columns allow null
            } else {
                _update(newValue);
            }
            function _addValidationCss() {
                // Show validation error
                if (settings.allowNulls.errorClass) {
                    $(inputField).addClass(settings.allowNulls.errorClass)
                } else {
                    $(inputField).css({ "border": "red solid 1px" })
                }
            }
            function _update(newValue) {
                var oldValue = cell.data();
                cell.data(newValue);
                //Return cell & row.
                settings.onUpdate(cell, row, oldValue);
            }

            //Redraw table
            table.draw();
        },
        // CANCEL
        cancelEditableCell: function (callingElement) {
            var cell = table.cell($(callingElement).parent());
            // Set cell to it's original value
            cell.data(cell.data())

            // Redraw table
            table.draw();
        }
    });



    // On cell click
    $(table.body()).on('click', 'td', function () {

        var currentColumnIndex = table.cell(this).index().column;

        // DETERMINE WHAT COLUMNS CAN BE EDITED
        if ((settings.columns && settings.columns.indexOf(currentColumnIndex) > -1) || (!settings.columns)) {
            var row = table.row($(this).parents('tr'));
            editableCellsRow = row;

            var cell = table.cell(this).node();

            // Show input
            if (!$(cell).find('input').length && !$(cell).find('select').length) {
                // Input CSS
                var input = getInputHtml(currentColumnIndex, settings);
                $(cell).html(input.html);
                if (input.focus) {
                    $('#ejbeatycelledit').focus();
                }
            }
        }
    });
});

function getInputHtml(currentColumnIndex, settings) {
    var inputSetting, inputType, input, inputCss, confirmCss, cancelCss;

    input = {"focus":true,"html":null}

    $.each(settings.inputTypes, function (index, setting) {
        if (setting.column == currentColumnIndex) {
            inputSetting = setting;
            inputType = inputSetting.type.toLowerCase();
        }
    });

    if (settings.inputCss) { inputCss = settings.inputCss; }
    if (settings.confirmationButton) {
        confirmCss = settings.confirmationButton.confirmCss;
        cancelCss = settings.confirmationButton.cancelCss;
        inputType = inputType + "-confirm";
    }
    switch (inputType) {
        case "list":
            input.html = "<select class='" + inputCss + "' onchange='$(this).updateEditableCell(this);'>";
            $.each(inputSetting.options, function (index, option) {
                input.html = input.html + "<option value='" + option.value + "' >" + option.display + "</option>"
            });
            input.html = input.html + "</select>";
            input.focus = false;
            break;
        case "list-confirm": // List w/ confirm
            input.html = "<select class='" + inputCss + "'>";
            $.each(inputSetting.options, function (index, option) {
                input.html = input.html + "<option value='" + option.value + "' >" + option.display + "</option>"
            });
            input.html = input.html + "</select>&nbsp;<a href='#' class='" + confirmCss + "' onclick='$(this).updateEditableCell(this);'>Confirm</a> <a href='#' class='" + cancelCss + "' onclick='$(this).cancelEditableCell(this)'>Cancel</a> ";
            input.focus = false;
            break;
        case "text-confirm": // text input w/ confirm
            input.html = "<input id='ejbeatycelledit' class='" + inputCss + "'>&nbsp;<a href='#' class='" + confirmCss + "' onclick='$(this).updateEditableCell(this)'>Confirm</a> <a href='#' class='" + cancelCss + "' onclick='$(this).cancelEditableCell(this)'>Cancel</a> ";
            break;
        case "undefined-confirm": // text input w/ confirm
            input.html = "<input id='ejbeatycelledit' class='" + inputCss + "'>&nbsp;<a href='#' class='" + confirmCss + "' onclick='$(this).updateEditableCell(this)'>Confirm</a> <a href='#' class='" + cancelCss + "' onclick='$(this).cancelEditableCell(this)'>Cancel</a> ";
            break;
        default: // text input
            input.html = "<input id='ejbeatycelledit' class='" + inputCss + "' onfocusout='$(this).updateEditableCell(this)'>";
            break;
    }
    return input;
}

function getInputField(callingElement) {
    // Update datatables cell value
    var inputField;
    switch ($(callingElement).prop('nodeName').toLowerCase()) {
        case 'a': // This means they're using confirmation buttons
            if ($(callingElement).siblings('input').length > 0) {
                inputField = $(callingElement).siblings('input');
            }
            if ($(callingElement).siblings('select').length > 0) {
                inputField = $(callingElement).siblings('select');
            }
            break;
        default:
            inputField = $(callingElement);
    }
    return inputField;
}
